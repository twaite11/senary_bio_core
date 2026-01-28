import subprocess
import os

class FusionHunter:
    def __init__(self, star_fusion_path="/usr/local/bin/STAR-Fusion"):
        self.tool_path = star_fusion_path
        self.genome_lib = "data/libraries/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play" # You download this once

    def run_analysis(self, left_fastq, right_fastq, output_dir):
        """
        Runs STAR-Fusion on raw paired-end FASTQ files to find novel chimeras.
        """
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        cmd = [
            self.tool_path,
            "--left_fq", left_fastq,
            "--right_fq", right_fastq,
            "--genome_lib_dir", self.genome_lib,
            "--output_dir", output_dir,
            "--CPU", "8" # Use all your cores
        ]

        print(f"[*] Running Fusion Detection on {left_fastq}...")
        try:
            subprocess.run(cmd, check=True)
            print(f"[SUCCESS] Fusion report generated in {output_dir}")
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] STAR-Fusion failed: {e}")

