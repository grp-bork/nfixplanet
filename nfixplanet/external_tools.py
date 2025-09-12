import shutil
import subprocess

def prodigal(input_file: str, out_file: str):
    """Run prodigal

    Arguments:
        input_file (str): full path of input fasta
        out_file (str): fullpath of output file
    """
    try:
        if input_file.endswith(".gz"):
            uncompressed = subprocess.Popen(
                ("zcat", input_file), stdout=subprocess.PIPE
            )
            subprocess.check_output(
                [
                    "prodigal",
                    "-i",
                    "/dev/stdin",
                    "-a",
                    out_file,
                    "-o",
                    "/dev/null",
                    "-p",
                    "meta",
                    "-q",
                ],
                stdin=uncompressed.stdout,
                universal_newlines=True,
            )
        else:
            subprocess.check_output(
                [
                    "prodigal",
                    "-i",
                    input_file,
                    "-a",
                    out_file,
                    "-o",
                    "/dev/null",
                    "-p",
                    "meta",
                    "-q",
                ],
                universal_newlines=True,
            )
    except subprocess.CalledProcessError:
        print(f"[ERROR] Failed to run Prodigal {input_file}")


def check_if_tool_exists(tool_name: str) -> bool:
    """Check if tool is available."""
    return shutil.which(tool_name) is not None
