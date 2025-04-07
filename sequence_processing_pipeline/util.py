import re


PAIR_UNDERSCORE = (re.compile(r'_R1_'), '_R1_', '_R2_')
PAIR_DOT = (re.compile(r'\.R1\.'), '.R1.', '.R2.')
PAIR_TESTS = (PAIR_UNDERSCORE, PAIR_DOT)


def determine_orientation(file_name):
    # aka forward, reverse, and indexed reads
    orientations = ['R1', 'R2', 'I1', 'I2']

    results = []

    # assume orientation is always present in the file's name.
    # assume that it is of one of the four forms above.
    # assume that it is always the right-most occurance of the four
    # orientations above.
    # assume that orientation is encapsulated with either '_' or '.'
    # e.g.: '_R1_', '.I2.'.
    # assume users can and will include any or all of the four
    # orientation as part of their filenames as well. e.g.:
    # ABC_7_04_1776_R1_SRE_S3_L007_R2_001.trimmed.fastq.gz
    for o in orientations:
        variations = [f"_{o}_", f".{o}."]
        for v in variations:
            # rfind searches from the end of the string, rather than
            # its beginning. It returns the position in the string
            # where the substring begins.
            results.append((file_name.rfind(v), o))

    # the orientation will be the substring found with the maximum
    # found value for pos. That is, it will be the substring that
    # begins at the rightest most position in the file name.
    results.sort(reverse=True)

    pos, orientation = results[0]

    # if no orientations were found, then return None.
    return None if pos == -1 else orientation


def iter_paired_files(files):
    """Yield matched r1/r2 paired files"""
    files = sorted(files)

    if len(files) % 2 != 0:
        raise ValueError("Files are not paired")

    files = iter(files)
    r1 = iter(files)
    r2 = iter(files)

    for r1_fp, r2_fp in zip(r1, r2):
        matched = False
        for pattern, r1_exp, r2_exp in PAIR_TESTS:
            if pattern.search(r1_fp):
                if r2_exp not in r2_fp:
                    raise ValueError(f"Cannot find '{r2_exp}' in '{r2_fp}'")

                # replace find w/find so that search for R1 and R2 begin
                # from the end of the string, not the beginning. This prevents
                # the code from breaking when filenames include R1 and R2 as
                # part of their name in addition to representing forward and
                # reversed reads e.g.:
                # LS_8_22_2014_R1_SRE_S3_L007_R1_001.trimmed.fastq.gz
                # LS_8_22_2014_R1_SRE_S3_L007_R2_001.trimmed.fastq.gz
                # using find(), r1_prefix and r2_prefix will be the following:
                # r1_prefix will be: LS_8_22_2014
                # r2_prefix will be: LS_8_22_2014_R1_SRE_S3_L007
                r1_prefix = r1_fp[:r1_fp.rfind(r1_exp)]
                r2_prefix = r2_fp[:r2_fp.rfind(r2_exp)]

                if r1_prefix != r2_prefix:
                    raise ValueError(f"Mismatch prefixes:\n{r1_prefix}\n"
                                     f"{r2_prefix}")

                matched = True
                break

        if not matched:
            raise ValueError(f"Unable to match:\n{r1_fp}\n{r2_fp}")

        yield (r1_fp, r2_fp)
