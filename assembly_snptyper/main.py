import argparse
import io
import logging
import multiprocessing
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

from assembly_snptyper.version import __version__


def check_external_dependencies():
    """
    Check if external dependencies are available in PATH
    """
    try:
        subprocess.run(
            ["minimap2", "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
    except FileNotFoundError:
        raise FileNotFoundError("minimap2 not found in PATH")

    try:
        subprocess.run(
            ["samtools", "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
    except FileNotFoundError:
        raise FileNotFoundError("samtools not found in PATH")


def convert_vcf_to_bed(vcf, bed_path):
    """
    Convert reference VCF to BED file

    The BED file is later used by samtools mpileup to filter variants

    Parameters
    ----------
    vcf : Path
        Path to reference VCF
    bed_path : Path
        Path to write BED file

    Returns
    -------
    None

    Raises
    ------
    AssertionError
        If REF or ALT field in VCF is not at most a single character long
    """
    # Expecting a single-sample VCF here
    data = pd.read_csv(
        vcf,
        comment="#",
        sep="\t",
        header=None,
        names=[
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            "SAMPLE",
        ],
        usecols=["#CHROM", "POS", "REF", "ALT"],
    )

    # Check if no indels are listed
    assert (
        data["REF"].str.len().max() == 1
    ), "REF field should be a single character: Reference VCF should not have indels"
    assert (
        data["ALT"].str.len().max() == 1
    ), "ALT field should be a single character: Reference VCF should not have indels"

    # Init Df which will be written to BED
    output = data[["#CHROM", "POS"]].copy()
    output.columns = ["CHROM", "END"]

    output["START"] = data["POS"].astype(int) - 1

    output[["CHROM", "START", "END"]].to_csv(
        bed_path, sep="\t", header=False, index=False
    )


def run_oneliner(bed_path, reference, input_asm, minimap_preset):
    """
    Run a minimap2 + samtools oneliner to generate a pileup from an input assembly

    Parameters
    ----------
    bed_path : Path
        Path to bed file
    reference : Path
        Path to reference genome
    input_asm : Path
        Path to single input assembly

    Returns
    -------
    str
        Output of samtools mpileup as text string
    """
    # run minimap2, pipe to samtools view, samtools sort and samtools mpileup
    # return stdout as text string
    minimap2 = subprocess.Popen(
        ["minimap2", "-ax", minimap_preset, reference, input_asm],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    samtools_view = subprocess.Popen(
        ["samtools", "view", "-b", "-"], stdin=minimap2.stdout, stdout=subprocess.PIPE
    )
    samtools_sort = subprocess.Popen(
        ["samtools", "sort", "-l", "0", "-"],
        stdin=samtools_view.stdout,
        stdout=subprocess.PIPE,
    )
    samtools_mpileup = subprocess.Popen(
        ["samtools", "mpileup", "-aa", "--positions", bed_path, "-f", reference, "-"],
        stdin=samtools_sort.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    output_text = samtools_mpileup.stdout.read().decode("utf-8")
    logging.debug(f"{output_text=}")
    return output_text


def parse_mpileup_output(mpileup_output, vcf, sample):
    """
    Parse the output of samtools mpileup and compare with the reference VCF.

    Prints the result to STDOUT

    Parameters
    ----------
    mpileup_output : str
        Output of samtools mpileup as text string
    vcf : Path
        Path to reference VCF
    sample : str
        Sample name

    Returns
    -------
    None

    Notes
    -----
    If a variant of the type of interest (ALT in reference VCF) has a depth >1,
    the variant might be located in a possibly duplicated region. The variant is
    still reported as matching in that case.

    """
    # Read data from text string
    data = pd.read_csv(
        io.StringIO(mpileup_output),
        header=None,
        sep="\t",
        names=["CHROM", "POS", "REF", "DEPTH", "ALT", "QUAL"],
    )

    # Expecting a single-sample VCF here
    ref = pd.read_csv(
        vcf,
        comment="#",
        sep="\t",
        header=None,
        names=[
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            "SAMPLE",
        ],
        usecols=["POS", "REF", "ALT"],
    )

    # Merge data
    merged = ref.merge(
        data, left_on="POS", right_on="POS", how="left", suffixes=("_ref", "_obs")
    )

    output = {}
    output["variants_missing"] = sum(merged["ALT_obs"] == "*")  # zero cov or deletion
    output["variants_multiple_cov"] = sum(
        merged["DEPTH"] > 1
    )  # cov >1: possible duplication
    output["variants_in_scheme"] = len(merged)  # total variants in reference vcf
    output["matching_variants"] = sum(
        merged["ALT_ref"].str.upper() == merged["ALT_obs"].str.upper()
    )  # variants matching type of interest
    output["wt_variants"] = sum(
        (merged["DEPTH"] == 1)
        & ((merged["ALT_obs"] == ".") | (merged["ALT_obs"] == ","))
    )  # valid (DP==1) variants that are wildtype (ALT==. or ,)

    return output


def wrapper(args_dict):
    """
    Wrapper function to run the typing workflow on a single input assembly

    This is needed to pass a single argument to the multiprocessing pool

    Parameters
    ----------
    args_dict : dict
        Dictionary with arguments for the workflow

    Returns
    -------
    None

    """
    sample = Path(args_dict["input_asm"]).stem
    result = run_oneliner(
        args_dict["bed_path"], args_dict["reference"], args_dict["input_asm"], args_dict["minimap_preset"]
    )
    output = parse_mpileup_output(result, args_dict["vcf"], sample)
    logging.info(f"Processed {sample}")
    output["sample"] = sample
    return output


def run_parallel(bed_path, reference, list_input, vcf, minimap_preset, num_processes):
    """
    Run the typing workflow in parallel on multiple input assemblies

    Parameters
    ----------
    bed_path : Path
        Path to bed file
    reference : Path
        Path to reference genome
    list_input : list
        List object of input assemblies
    vcf : Path
        Path to reference VCF
    num_processes : int
        Number of processes to run in parallel

    Returns
    -------
    None

    """
    list_wrapper_args = [
        {
            "bed_path": bed_path,
            "reference": reference,
            "input_asm": input_asm,
            "vcf": vcf,
            "minimap_preset": minimap_preset,
        }
        for input_asm in list_input
    ]
    logging.info(
        f"Running {len(list_wrapper_args)} samples in parallel with {num_processes} processes"
    )
    logging.info("Starting multiprocessing")
    pool = multiprocessing.Pool(num_processes)
    list_of_outputs = pool.map(wrapper, list_wrapper_args)
    pool.close()
    pool.join()
    logging.info("Finished multiprocessing")
    logging.info("Printing results")
    list_cols = [
        "sample",
        "matching_variants",
        "wt_variants",
        "variants_in_scheme",
        "variants_missing",
        "variants_multiple_cov",
    ]
    pd.DataFrame(list_of_outputs)[list_cols].to_csv(sys.stdout, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--vcf",
        required=True,
        type=Path,
        help="VCF file with variants that determine the type",
    )
    parser.add_argument(
        "--reference", required=True, type=Path, help="Reference genome"
    )
    parser.add_argument(
        "--list_input", required=True, type=Path, help="List of input assemblies"
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        help="Number of processes passed to multiprocessing",
        default=1,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Verbose output, can be used multiple times. 0 = warning, 1 = info, 2 = debug",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    args = parser.parse_args()

    # set up logger
    if args.verbose == 0:
        log_level = logging.WARNING
    elif args.verbose == 1:
        log_level = logging.INFO
    elif args.verbose >= 2:
        log_level = logging.DEBUG

    logging.basicConfig(
        level=log_level,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="[%Y-%m-%d %H:%M:%S]",
    )

    # Not exposed in CLI yet, needs testing
    minimap_preset = "asm5"
    logging.info(f"Using minimap2 preset: {minimap_preset}")

    # check non-python dependencies
    logging.info("Checking external dependencies")
    check_external_dependencies()
    logging.info("External dependencies found")

    # read list of input assemblies
    logging.info(f"Reading list of input assemblies from {args.list_input}")
    with open(args.list_input) as f:
        list_input = [Path(line.strip()) for line in f]
    logging.info(f"Number of input assemblies: {len(list_input)}")

    # create temporary bed file, is not deleted on closing of script
    # this is possible in python3.12 with delete_on_close=False + delete=True
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".bed") as tmp:
        bed_path = Path(tmp.name)

    # make bed which is used by samtools mpileup to filter
    convert_vcf_to_bed(args.vcf, bed_path)
    logging.info(f"Created temporary bed file: {bed_path}")

    logging.info("Starting typing workflow")
    run_parallel(bed_path, args.reference, list_input, args.vcf, minimap_preset, args.processes)


if __name__ == "__main__":
    main()
