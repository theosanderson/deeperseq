# use argparse to get a single filename
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("filename")
args = parser.parse_args()
filename = args.filename

# read from stdin, write to stdout, count lines and write the number to filename every 10000 lines
import sys
count = 0
for line in sys.stdin:
    count += 1
    if count % 10000 == 0:
        with open(filename, "w") as f:
            f.write(str(count))
    sys.stdout.write(line)