#!/usr/bin/env python3
"""
Generates the genome-wide aneuploidy ideogram, circular plot and text
summary.

Source: {sample}_cnv_bands.tsv (export_cnv_tables.py) -- bands already
classified (gain/loss/neutral) by sex-aware absolute cn, not an
independent log2 cutoff. Same source/classification used by the scatter
plots, so color and findings agree across plot types.

Usage:
    python plot_aneuploidy_ideogram.py \
        --bands karyogram/sample/sample_cnv_bands.tsv \
        --sample ERR2990113 \
        --output karyogram/sample/sample_aneuploidy_ideogram.pdf
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import numpy as np
import argparse
import logging
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cnvkaryo_plot_common import CHR_ORDER, CHR_LENGTHS, GAIN_COLOR, LOSS_COLOR

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)

COLOR_GAIN = GAIN_COLOR
COLOR_LOSS = LOSS_COLOR
COLOR_NEUTRAL = '#ecf0f1'

# Vertical split of each chromosome row: call on top, band detail below.
BAND_HEIGHT = 0.22
BAND_TOP = 0.30


class AneuploidyPlotter:
    def __init__(self, sample_id):
        self.sample_id = sample_id
        self.chr_order = CHR_ORDER

    def load_cnvs(self, bands_file):
        """Loads already-classified bands (gain/loss/neutral) from
        {sample}_cnv_bands.tsv -- does not reclassify by log2 here; uses
        the classification already computed (absolute cn vs sex-aware
        expected ploidy) in export_cnv_tables.py."""
        logger.info(f"Loading bands: {bands_file}")

        bands = pd.read_csv(bands_file, sep='\t', comment='#')
        bands = bands.rename(columns={
            'classification': 'cnv_type',
            'log2_weighted_mean': 'log2',
        })

        logger.info(f"  Total: {len(bands)}, Altered: {(bands['cnv_type'] != 'neutral').sum()}")

        return bands

    def load_arms(self, arms_file):
        """Whole-chromosome and arm calls from call_arm_level.py. These, not
        the CBS bands, decide what the figure reports: segmentation
        optimises breakpoint placement, so on noisy data it produces many
        short high-variance segments and a chromosome-scale verdict read
        off their total length is wrong. Measured on this panel, the band
        path reported 156 non-neutral bands across 17 chromosomes for a
        constitutional trisomy 21, whose only true abnormality is +21."""
        logger.info(f"Loading arm-level calls: {arms_file}")
        # iscn as str: pandas reads "+21" as the float 21.0 otherwise, and
        # the karyotype label silently stops being a karyotype label.
        arms = pd.read_csv(arms_file, sep='\t', dtype={'iscn': str, 'region': str},
                           keep_default_na=False, na_values=[])
        called = arms[arms['call'].isin(['gain', 'loss'])]
        logger.info(f"  Regions: {len(arms)}, called: {len(called)}")
        return arms

    def called_regions(self, arms, chrom):
        """(start, end, call) for every called region on this chromosome. A
        whole-chromosome call subsumes its arms upstream, so these never
        overlap."""
        sel = arms[(arms['chrom'] == chrom) & (arms['call'].isin(['gain', 'loss']))]
        return [(int(r['start']), int(r['end']), r['call']) for _, r in sel.iterrows()]

    def calculate_chromosome_aneuploidy(self, arms):
        """Status per chromosome, taken straight from the caller's verdict
        rather than re-derived from how much of the chromosome is covered
        by altered segments."""
        chr_status = {}
        for chrom in self.chr_order:
            if chrom not in CHR_LENGTHS:
                continue
            sel = arms[arms['chrom'] == chrom]
            whole = sel[(sel['region'] == 'whole') & (sel['call'].isin(['gain', 'loss']))]
            if len(whole):
                call = whole.iloc[0]['call']
                chr_status[chrom] = 'trisomy' if call == 'gain' else 'monosomy'
                continue
            partial = sel[(sel['region'] != 'whole') & (sel['call'].isin(['gain', 'loss']))]
            if len(partial):
                calls = set(partial['call'])
                chr_status[chrom] = ('partial_gain' if calls == {'gain'}
                                     else 'partial_loss' if calls == {'loss'}
                                     else 'partial_mixed')
            else:
                chr_status[chrom] = 'neutral'
        return chr_status

    def plot_genome_ideogram(self, cnvs, arms, output_file=None):
        logger.info("Generating ideogram...")

        fig, axes = plt.subplots(
            nrows=24, ncols=1,
            figsize=(16, 20),
            gridspec_kw={'hspace': 0.3}
        )

        chr_status = self.calculate_chromosome_aneuploidy(arms)

        for idx, chrom in enumerate(self.chr_order):
            ax = axes[idx]

            if chrom not in CHR_LENGTHS:
                ax.axis('off')
                continue

            chr_length = CHR_LENGTHS[chrom]
            chr_cnvs = cnvs[cnvs['chrom'] == chrom]

            # Upper body: the karyotype call. Lower strip: the CBS bands
            # behind it, kept visible so a reader can see the underlying
            # evidence without it driving the verdict.
            ax.add_patch(Rectangle(
                (0, BAND_TOP), chr_length / 1e6, 1 - BAND_TOP,
                facecolor=COLOR_NEUTRAL, edgecolor='black', linewidth=0.5
            ))
            for start, end, call in self.called_regions(arms, chrom):
                ax.add_patch(Rectangle(
                    (start / 1e6, BAND_TOP), (end - start) / 1e6, 1 - BAND_TOP,
                    facecolor=COLOR_GAIN if call == 'gain' else COLOR_LOSS,
                    alpha=0.85, edgecolor='none'
                ))

            ax.add_patch(Rectangle(
                (0, 0), chr_length / 1e6, BAND_HEIGHT,
                facecolor='#ffffff', edgecolor='#999999', linewidth=0.4
            ))
            # Bands are coloured only inside a called region. Elsewhere they
            # are grey: the strip shows the evidence behind the call without
            # asserting alterations the caller rejected.
            for _, cnv in chr_cnvs.iterrows():
                if cnv['cnv_type'] not in ('gain', 'loss'):
                    continue
                ax.add_patch(Rectangle(
                    (cnv['start'] / 1e6, 0), (cnv['end'] - cnv['start']) / 1e6, BAND_HEIGHT,
                    facecolor='#999999', alpha=0.35, edgecolor='none'
                ))

            # Common scale across all rows: drawing each chromosome to the
            # full panel width would show chr21 as long as chr1, which is
            # the one thing a karyogram must get right.
            ax.set_xlim(0, max(CHR_LENGTHS.values()) / 1e6)
            ax.set_ylim(0, 1)
            ax.set_yticks([])

            status = chr_status.get(chrom, 'neutral')
            status_label = ''

            iscn = ', '.join(
                arms[(arms['chrom'] == chrom)
                     & (arms['iscn'].astype(str) != '')]['iscn'].astype(str))
            if status == 'trisomy':
                status_label = f' [{iscn}]'
                ax.set_facecolor('#ffe6e6')
            elif status == 'monosomy':
                status_label = f' [{iscn}]'
                ax.set_facecolor('#e6f2ff')
            elif status.startswith('partial'):
                status_label = f' [{iscn}]'

            chr_label = chrom.replace('chr', '')
            ax.set_ylabel(f'{chr_label}{status_label}', fontsize=10, fontweight='bold')

            if idx == len(self.chr_order) - 1:
                ax.set_xlabel('Position (Mb)', fontsize=12)
            else:
                ax.set_xticklabels([])

        fig.suptitle(
            f'Aneuploidy Ideogram - {self.sample_id}',
            fontsize=16, fontweight='bold', y=0.995
        )

        legend_elements = [
            mpatches.Patch(facecolor=COLOR_GAIN, label='Gain'),
            mpatches.Patch(facecolor=COLOR_LOSS, label='Loss'),
            mpatches.Patch(facecolor='#999999',
                           label='grey: CBS evidence the caller did not call'),
        ]

        fig.legend(
            handles=legend_elements,
            loc='upper right',
            bbox_to_anchor=(0.98, 0.99),
            fontsize=10
        )

        plt.tight_layout()

        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"Ideogram saved: {output_file}")

        return fig

    def plot_cnv_summary_circular(self, cnvs, output_file=None):
        logger.info("Generating circular plot...")

        fig = plt.figure(figsize=(14, 14))
        ax = fig.add_subplot(111, projection='polar')

        total_length = sum(CHR_LENGTHS.values())
        chr_angles = {}
        cumulative = 0

        for chrom in self.chr_order:
            if chrom not in CHR_LENGTHS:
                continue

            start_angle = (cumulative / total_length) * 2 * np.pi
            chr_length = CHR_LENGTHS[chrom]
            end_angle = ((cumulative + chr_length) / total_length) * 2 * np.pi

            chr_angles[chrom] = {
                'start': start_angle,
                'end': end_angle,
                'mid': (start_angle + end_angle) / 2
            }

            cumulative += chr_length

        for chrom, angles in chr_angles.items():
            theta = np.linspace(angles['start'], angles['end'], 100)
            r = np.ones_like(theta) * 1.0

            ax.plot(theta, r, color='black', linewidth=2)
            ax.fill_between(theta, 0.95, 1.05, alpha=0.1, color='gray')

            chr_label = chrom.replace('chr', '')
            ax.text(
                angles['mid'], 1.15, chr_label,
                ha='center', va='center', fontsize=10, fontweight='bold'
            )

        for _, cnv in cnvs.iterrows():
            chrom = cnv['chrom']

            if chrom not in chr_angles:
                continue

            chr_info = chr_angles[chrom]
            chr_length = CHR_LENGTHS[chrom]

            start_frac = cnv['start'] / chr_length
            end_frac = cnv['end'] / chr_length

            theta_start = chr_info['start'] + start_frac * (chr_info['end'] - chr_info['start'])
            theta_end = chr_info['start'] + end_frac * (chr_info['end'] - chr_info['start'])

            theta_cnv = np.linspace(theta_start, theta_end, 50)

            if cnv['cnv_type'] == 'gain':
                r_cnv = np.ones_like(theta_cnv) * (1.0 + cnv['log2'] * 0.3)
                color = COLOR_GAIN
            elif cnv['cnv_type'] == 'loss':
                r_cnv = np.ones_like(theta_cnv) * (1.0 + cnv['log2'] * 0.3)
                color = COLOR_LOSS
            else:
                continue

            ax.plot(theta_cnv, r_cnv, color=color, linewidth=2, alpha=0.7)

        ax.set_ylim(0, 1.5)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.spines['polar'].set_visible(False)
        ax.grid(False)

        plt.title(
            f'Circular CNV Profile - {self.sample_id}',
            fontsize=16, fontweight='bold', pad=20
        )

        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"Circular plot saved: {output_file}")

        return fig

    def generate_aneuploidy_summary(self, cnvs, arms):
        chr_status = self.calculate_chromosome_aneuploidy(arms)

        summary = {
            'trisomies': [],
            'monosomies': [],
            'partial_gains': [],
            'partial_losses': [],
            'aneuploidy_score': 0
        }

        for chrom, status in chr_status.items():
            if status == 'trisomy':
                summary['trisomies'].append(chrom)
            elif status == 'monosomy':
                summary['monosomies'].append(chrom)
            elif status == 'partial_gain':
                summary['partial_gains'].append(chrom)
            elif status == 'partial_loss':
                summary['partial_losses'].append(chrom)

        # Scored over called regions, so the figure's headline number and
        # its colours come from the same source.
        called = arms[arms['call'].isin(['gain', 'loss'])]
        altered_bases = int((called['end'] - called['start']).sum())
        total_bases = sum(CHR_LENGTHS[c] for c in self.chr_order
                          if c in CHR_LENGTHS and c not in ('chrX', 'chrY'))
        summary['aneuploidy_score'] = altered_bases / total_bases if total_bases else 0
        summary['findings'] = [x for x in called['iscn'].astype(str).tolist() if x]
        return summary


def main():
    parser = argparse.ArgumentParser(
        description='Generates the genome-wide aneuploidy ideogram'
    )
    parser.add_argument('--bands', required=True, help='{sample}_cnv_bands.tsv (export_cnv_tables.py) -- detail layer')
    parser.add_argument('--arms', required=True, help='{sample}_arm_level.tsv (call_arm_level.py) -- the karyotype call')
    parser.add_argument('--sample', required=True, help='Sample ID')
    parser.add_argument('--output', required=True, help='Output PDF file')
    parser.add_argument('--output-circular', help='Circular plot (optional)')
    parser.add_argument('--output-summary', help='Text summary (optional)')

    args = parser.parse_args()

    plotter = AneuploidyPlotter(args.sample)
    cnvs = plotter.load_cnvs(args.bands)
    arms = plotter.load_arms(args.arms)
    plotter.plot_genome_ideogram(cnvs, arms, args.output)

    if args.output_circular:
        plotter.plot_cnv_summary_circular(cnvs, args.output_circular)

    summary = plotter.generate_aneuploidy_summary(cnvs, arms)

    logger.info("\n" + "=" * 60)
    logger.info("ANEUPLOIDY SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Aneuploidy score: {summary['aneuploidy_score']:.2%}")

    if summary['trisomies']:
        logger.info(f"\nTrisomies: {', '.join(summary['trisomies'])}")
    if summary['monosomies']:
        logger.info(f"Monosomies: {', '.join(summary['monosomies'])}")
    if summary['partial_gains']:
        logger.info(f"Partial gains: {', '.join(summary['partial_gains'])}")
    if summary['partial_losses']:
        logger.info(f"Partial losses: {', '.join(summary['partial_losses'])}")

    if args.output_summary:
        with open(args.output_summary, 'w') as f:
            f.write(f"Sample: {args.sample}\n")
            f.write(f"Findings (ISCN): "
                    f"{', '.join(summary['findings']) if summary['findings'] else 'None'}\n")
            f.write(f"Fraction of autosomal sequence called: "
                    f"{summary['aneuploidy_score']:.4f}\n\n")
            f.write(f"Trisomies: {', '.join(summary['trisomies']) if summary['trisomies'] else 'None'}\n")
            f.write(f"Monosomies: {', '.join(summary['monosomies']) if summary['monosomies'] else 'None'}\n")
            f.write(f"Partial Gains: {', '.join(summary['partial_gains']) if summary['partial_gains'] else 'None'}\n")
            f.write(f"Partial Losses: {', '.join(summary['partial_losses']) if summary['partial_losses'] else 'None'}\n")

        logger.info(f"\nSummary saved: {args.output_summary}")

    logger.info("\nPlots generated successfully.")


if __name__ == '__main__':
    main()
