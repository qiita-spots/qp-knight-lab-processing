from os import walk
from sys import argv
from os.path import join, split
from json import dumps


def extract_metadata(log_output_file_path):
    with open(log_output_file_path, 'r') as f:
        lines = f.readlines()
        lines = [x.strip() for x in lines]
        if len(lines) != 2:
            raise ValueError("error processing %s" % log_output_file_path)
        _dir, _file = split(lines[0])
        seq_counts, base_pairs = lines[1].split('\t')
        return _dir, _file, int(seq_counts), int(base_pairs)


def aggregate_counts(fp):
    results = {}

    for root, dirs, files in walk(fp):
        for _file in files:
            if _file.endswith('.out'):
                log_output_file = join(root, _file)
                _dir, _file, seq_counts, base_pairs = \
                    extract_metadata(log_output_file)

                if _dir not in results:
                    results[_dir] = {}

                results[_dir][_file] = {'seq_counts': seq_counts,
                                        'base_pairs': base_pairs}

    return results


if __name__ == '__main__':
    results = aggregate_counts(argv[1])
    with open(argv[2], 'w') as f:
        print(dumps(results, indent=2), file=f)
