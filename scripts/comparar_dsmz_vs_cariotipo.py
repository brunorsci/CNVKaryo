#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Confronta o cariotipo em bandas do DSMZ com o resultado de profundidade.

CRITERIO DE PONTUACAO, FIXADO ANTES DE MEDIR
============================================
Um cariotipo em bandas NAO e padrao-ouro irrestrito para um metodo baseado em
profundidade de leitura, e usa-lo como tal produz numero bonito e errado. Tres
razoes, e cada uma vira uma regra:

1. BANDAS NAO VEEM DOSAGEM EM TRANSLOCACAO EQUILIBRADA. t(1;19) nao muda
   numero de copias e nao pode ser exigida de um metodo de cobertura.
   -> So entram eventos que MUDAM DOSAGEM.

2. RESOLUCAO. Bandeamento resolve ~5-10 Mb; a analise por braco aqui usa a
   MEDIANA do braco e exige 50 janelas. Um del(7)(p21.1p21.2), que tem alguns
   Mb, nao desloca a mediana de 7p e nao deveria contar como falha.
   -> So entram eventos que cobrem >= 50% do braco, medido nas coordenadas das
      bandas. Os demais ficam numa categoria propria, "abaixo da escala do
      metodo", relatada a parte e nao somada como erro.

3. AUSENCIA NO CARIOTIPO NAO E AUSENCIA DE CNV. O bandeamento nao enxerga o
   que esta abaixo da sua resolucao, e o estoque do DSMZ nao e o mesmo estoque
   que foi sequenciado -- linhagem deriva com a passagem.
   -> Calcula-se RECALL (sensibilidade) contra os eventos documentados.
      NAO se calcula especificidade nem VPP contra o DSMZ, porque isso exigiria
      tratar o silencio do cariotipo como negativo verdadeiro, o que ele nao e.

Eventos ambiguos quanto a dosagem -- der(), add(), mar, t() -- ficam fora dos
dois lados da conta: nao sao exigidos nem contados contra.

4. CONSEQUENCIA DA REGRA 3, QUE E FACIL DEIXAR PASSAR: em cariotipo complexo,
   "-N" raramente significa que o cromossomo N desapareceu. O material costuma
   estar retido num derivativo, que o MESMO cariotipo lista a parte, e entao so
   um dos bracos e de fato perdido. Medido em 28/10/2026: BV-173 tem "-9" com
   9p em -0,006 e 9q em -0,74; HL-60 e GRANTA-519 tem "-17" com 17p em -0,67 e
   17q neutro; GRANTA-519 tem "-18" com 18p em -0,92 e 18q em +0,28. Exigir os
   dois bracos nesses casos nao mede o metodo, mede a convencao da notacao.
   -> "+N" e "-N" so entram quando o cariotipo NAO traz der()/t()/add()/i()/mar
      envolvendo aquele mesmo cromossomo. Havendo, a dosagem por braco e
      indeterminada e o evento sai da conta, como qualquer outro ambiguo.

Entrada: a tabela do DSMZ, o mapeamento amostra->linhagem, a tabela de bracos
do workflow e o arquivo de bandas. Saida: uma linha por evento confrontado.
"""
import argparse
import csv
import re
import sys
from collections import defaultdict

# Eventos de dosagem nao ambiguos. Deliberadamente nao se tenta inferir
# dosagem de der()/add()/mar: o desequilibrio depende do parceiro, que a
# notacao nem sempre informa.
RE_GAIN_WHOLE = re.compile(r'(?<![\w(])\+(\d{1,2})(?![\w)\d])')
RE_LOSS_WHOLE = re.compile(r'(?<![\w(])-(\d{1,2})(?![\w)\d])')
RE_DEL = re.compile(r'\bdel\((\d{1,2}|[XY])\)\(([pq][0-9.]*)(?:([pq][0-9.]*))?\)', re.I)
RE_DUP = re.compile(r'\bdup\((\d{1,2}|[XY])\)\(([pq][0-9.]*)(?:([pq][0-9.]*))?\)', re.I)
RE_ISO = re.compile(r'\bi\((\d{1,2}|[XY])\)\(([pq])10\)', re.I)


def carregar_bandas(caminho):
    bandas = defaultdict(list)
    with open(caminho) as fh:
        primeira = fh.readline()
        if not primeira.lower().startswith(('contig', 'chrom', '#')):
            fh.seek(0)
        for linha in fh:
            c = linha.rstrip('\n').split('\t')
            if len(c) < 5:
                continue
            chrom = c[0] if c[0].startswith('chr') else f'chr{c[0]}'
            try:
                bandas[chrom].append((int(c[1]), int(c[2]), c[3], c[4]))
            except ValueError:
                continue
    for c in bandas:
        bandas[c].sort()
    return bandas


def limites_braco(bandas, chrom):
    bs = bandas.get(chrom)
    if not bs:
        return None
    acen = [b for b in bs if b[3] == 'acen']
    if not acen:
        return None
    ini, fim = min(b[0] for b in bs), max(b[1] for b in bs)
    return {'p': (ini, min(b[0] for b in acen)),
            'q': (max(b[1] for b in acen), fim)}


def posicao_banda(bandas, chrom, nome):
    """(inicio, fim) da banda; aceita prefixo (q21 casa com q21.1 se q21 exato
    nao existir, como e habitual na notacao resumida do DSMZ)."""
    bs = bandas.get(chrom, [])
    exatas = [b for b in bs if b[2] == nome]
    if exatas:
        return (min(b[0] for b in exatas), max(b[1] for b in exatas))
    pref = [b for b in bs if b[2].startswith(nome)]
    if pref:
        return (min(b[0] for b in pref), max(b[1] for b in pref))
    return None


def extensao(bandas, chrom, b1, b2):
    """Trecho afetado por del/dup. Com UMA banda a notacao significa evento
    terminal (daquela banda ate o telomero do braco); com DUAS, intersticial."""
    arm = b1[0]
    lim = limites_braco(bandas, chrom)
    if not lim or arm not in lim:
        return None
    p1 = posicao_banda(bandas, chrom, b1)
    if not p1:
        return None
    if b2:
        p2 = posicao_banda(bandas, chrom, b2)
        if not p2:
            return None
        return (min(p1[0], p2[0]), max(p1[1], p2[1]))
    ini_br, fim_br = lim[arm]
    return (ini_br, p1[1]) if arm == 'p' else (p1[0], fim_br)


RE_REARRANJO = re.compile(r'\b(?:der|t|add|i|ins|dic|idic|r|inv)\(([^)]*)\)', re.I)


def cromossomos_rearranjados(karyotype):
    """Cromossomos citados em qualquer rearranjo estrutural. Para estes, um
    "+N"/"-N" no mesmo cariotipo nao determina dosagem por braco -- ver a
    regra 4 do cabecalho."""
    envolvidos = set()
    for m in RE_REARRANJO.finditer(karyotype):
        for n in re.findall(r'\d{1,2}|[XY]', m.group(1)):
            envolvidos.add(f'chr{n}')
    if re.search(r'\bmar\b', karyotype, re.I):
        envolvidos.add('*')          # marcador de origem nao declarada
    return envolvidos


def eventos_dsmz(karyotype, bandas):
    """[(chrom, braco_ou_None, 'ganho'/'perda', fracao_do_braco, notacao)]"""
    # Corta o texto livre que o DSMZ anexa depois do cariotipo.
    for marca in ('Molec. Genetics', 'carries ', 'resembles ', 'related to ',
                  'expression of', 'confirmed by'):
        i = karyotype.find(marca)
        if i > 0:
            karyotype = karyotype[:i]
    ev = []
    rearr = cromossomos_rearranjados(karyotype)
    for regex, tipo in ((RE_GAIN_WHOLE, 'ganho'), (RE_LOSS_WHOLE, 'perda')):
        for m in regex.finditer(karyotype):
            chrom = f'chr{m.group(1)}'
            if chrom in rearr:
                continue          # dosagem por braco indeterminada (regra 4)
            ev.append((chrom, None, tipo, 1.0, m.group(0)))
    for regex, tipo in ((RE_DEL, 'perda'), (RE_DUP, 'ganho')):
        for m in regex.finditer(karyotype):
            chrom = f'chr{m.group(1)}'
            tr = extensao(bandas, chrom, m.group(2), m.group(3))
            lim = limites_braco(bandas, chrom)
            if not tr or not lim:
                continue
            arm = m.group(2)[0]
            ini, fim = lim[arm]
            frac = (tr[1] - tr[0]) / (fim - ini) if fim > ini else 0.0
            ev.append((chrom, arm, tipo, min(frac, 1.0), m.group(0)))
    for m in RE_ISO.finditer(karyotype):
        chrom = f'chr{m.group(1)}'
        mantido = m.group(2)
        ev.append((chrom, 'q' if mantido == 'p' else 'p', 'perda', 1.0, m.group(0)))
        ev.append((chrom, mantido, 'ganho', 1.0, m.group(0)))
    # O texto do DSMZ as vezes repete o mesmo evento (sublinhas, "x2"); cada
    # evento distinto conta uma vez.
    vistos, unicos = set(), []
    for e in ev:
        chave = (e[0], e[1], e[2])
        if chave not in vistos:
            vistos.add(chave)
            unicos.append(e)
    return unicos


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dsmz', required=True)
    ap.add_argument('--mapeamento', required=True)
    ap.add_argument('--bracos', required=True, help='all_samples_arm_level.tsv')
    ap.add_argument('--cytobands', required=True)
    ap.add_argument('--fracao-minima', type=float, default=0.50,
                    help='fracao do braco que o evento precisa cobrir para ser '
                         'exigivel da analise por braco (ver cabecalho)')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    bandas = carregar_bandas(args.cytobands)
    linhagem = {}
    for r in csv.DictReader(open(args.mapeamento), delimiter='\t'):
        linhagem[r['Cell_line']] = r['DNA_sample']
    dsmz = {r['cell_line']: r for r in csv.DictReader(open(args.dsmz), delimiter='\t')}

    # Chamadas do workflow, por (amostra, chrom, braco)
    chamado = {}
    for r in csv.DictReader(open(args.bracos), delimiter='\t'):
        cl = r['classe_do_braco']
        chamado[(r['sample'], r['chrom'], r['braco'][-1])] = (
            'ganho' if cl == 'ganho_de_braco' else
            'perda' if cl == 'perda_de_braco' else 'neutro',
            float(r['desvio_ajustado']), float(r['fracao_janelas_alteradas']))

    linhas = []
    for cell, r in sorted(dsmz.items()):
        amostra = linhagem.get(cell)
        if not amostra:
            continue
        for chrom, arm, tipo, frac, nota in eventos_dsmz(r['karyotype'], bandas):
            lim = limites_braco(bandas, chrom)
            if not lim:
                continue
            bracos = [arm] if arm else [a for a in ('p', 'q') if a in lim]
            for a in bracos:
                k = (amostra, chrom, a)
                got, desvio, fr = chamado.get(k, ('sem_dado', float('nan'), float('nan')))
                escala = 'sim' if frac >= args.fracao_minima else 'nao'
                if escala == 'nao':
                    res = 'abaixo_da_escala'
                elif got == tipo:
                    res = 'concordante'
                elif got == 'neutro':
                    res = 'nao_detectado'
                elif got == 'sem_dado':
                    res = 'sem_dado'
                else:
                    res = 'sentido_oposto'
                linhas.append(dict(
                    linhagem=cell, amostra=amostra, notacao_dsmz=nota,
                    chrom=chrom, braco=a, esperado=tipo,
                    fracao_do_braco=round(frac, 2), na_escala_do_metodo=escala,
                    obtido=got, desvio_ajustado=round(desvio, 3) if desvio == desvio else '',
                    fracao_janelas=round(fr, 2) if fr == fr else '',
                    resultado=res,
                    ploidia=r['ploidy'], numero_modal=r['modal_number']))

    campos = ['linhagem', 'amostra', 'notacao_dsmz', 'chrom', 'braco', 'esperado',
              'fracao_do_braco', 'na_escala_do_metodo', 'obtido', 'desvio_ajustado',
              'fracao_janelas', 'resultado', 'ploidia', 'numero_modal']
    with open(args.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=campos, delimiter='\t')
        w.writeheader()
        w.writerows(linhas)

    naescala = [l for l in linhas if l['na_escala_do_metodo'] == 'sim']
    conc = sum(1 for l in naescala if l['resultado'] == 'concordante')
    print(f"eventos de dosagem extraidos do DSMZ : {len(linhas)}")
    print(f"  na escala do metodo (>= {args.fracao_minima:.0%} do braco): {len(naescala)}")
    print(f"  abaixo da escala, relatados a parte  : {len(linhas) - len(naescala)}")
    if naescala:
        print(f"\nRECALL na escala do metodo: {conc}/{len(naescala)} = "
              f"{conc / len(naescala):.3f}")
        for res in ('concordante', 'nao_detectado', 'sentido_oposto', 'sem_dado'):
            n = sum(1 for l in naescala if l['resultado'] == res)
            if n:
                print(f"  {res:18} {n:4}")
    print(f"\n-> {args.out}")


if __name__ == '__main__':
    main()
