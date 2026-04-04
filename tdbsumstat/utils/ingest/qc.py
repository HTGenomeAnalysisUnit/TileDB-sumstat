"""GWAS-lab quality-control step (optional)."""
import logging
import os

import polars as pl

logger = logging.getLogger(__name__)


class QCMixin:
    """Mixin that provides optional gwaslab-based QC on ``self.chunk_pl``."""

    def qc_sumstat(self, file_path: str) -> None:
        """Run gwaslab QC checks on the harmonised data in ``self.chunk_pl``.

        Logs are saved to ``{self.uri}_logs/{file_stem}.log``.
        On success, ``self.chunk_pl`` is replaced with the QC-filtered
        gwaslab DataFrame (converted back to Polars).

        Parameters
        ----------
        file_path:
            Original file path; used only to derive a log file name.
        """
        import gwaslab as gl

        directory = self.uri + "_logs"
        filename = os.path.basename(file_path)
        file_name = filename.split(".")[0]

        os.makedirs(directory, exist_ok=True)

        sumstat_preqc = self.chunk_pl.to_pandas()

        if self.type_sumstat == "gwas":
            if self.type_trait == "quant":
                sumstat_gl = gl.Sumstats(
                    sumstat_preqc,
                    snpid="SNPID", chrom="CHR", pos="POS",
                    eaf="EAF", beta="BETA", se="SE", p="P", n="N",
                    ea="EA", nea="NEA",
                    other=["TRAIT", "RSID"],
                )
            else:
                sumstat_gl = gl.Sumstats(
                    sumstat_preqc,
                    snpid="SNPID", chrom="CHR", pos="POS",
                    eaf="EAF", beta="BETA", se="SE", p="P", n="N",
                    ncase="N_CASES", ncontrol="N_CONTROLS",
                    ea="EA", nea="NEA",
                    other=["TRAIT", "RSID"],
                )
        else:
            sumstat_gl = gl.Sumstats(
                sumstat_preqc,
                snpid="SNPID", chrom="CHR", pos="POS",
                eaf="EAF", beta="BETA", se="SE", p="P", n="N",
                other=["CELL", "GENE", "RSID", "DIST", "PHENO_VAR"],
            )

        sumstat_gl.fix_chr(remove=True)
        sumstat_gl.fix_pos(remove=True)
        sumstat_gl.fix_allele(remove=True)
        sumstat_gl.check_sanity()
        sumstat_gl.check_data_consistency()

        sumstat_gl.log.save(os.path.join(directory, file_name))
        self.chunk_pl = pl.from_pandas(sumstat_gl.data)
        logger.info("QC completed for %s", file_path)
