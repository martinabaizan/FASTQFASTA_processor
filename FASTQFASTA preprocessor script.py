import sys
import argparse
from sys import exit
import numbers
import math

def stats_update(stats, trimmed_stats, sequence, trimmed_seq):
    #Stats for bases processed
    stats["reads"] += 1
    stats["bases"] += len(sequence)
    for nt in sequence:
        stats[nt] += 1

    #Stats for bases trimmed
    if args.trim_left and args.trim_right is not None:
         trimmed_stats["bases"] += args.trim_left + args.trim_right
    elif args.trim_left is not None:
        trimmed_stats["bases"] += args.trim_left
    elif args.trim_right is not None:
        trimmed_stats["bases"] += args.trim_right

    if args.trim_left and args.trim_right is not None:
         trimmed_seq = (sequence[: int(args.trim_left)]) + (sequence[((len(sequence)) - args.trim_right):])
    elif args.trim_left is not None:
        trimmed_seq = sequence[: int(args.trim_left)]
    elif args.trim_right is not None:
        trimmed_seq = sequence[((len(sequence)) - args.trim_right):]
    for nt in trimmed_seq:
        trimmed_stats[nt] += 1

def print_summary(stats, trimmed_stats=None):
    # Translate operation
    if args.operation == "rc":
        operation = "reversed-complemented"
    elif args.operation == "trim":
        operation = "hard-trimmed"
    elif args.operation == "adaptor-removal":
        operation = "processed"
    # Print summary
    print("File '%s' has been successfully %s ('%s')" % (args.input, operation, args.output))
    print("Summary:")
    # total number of reads processed (i.e sequences)
    print("\t%d reads processed" % stats["reads"])
    # total number of bases processed (together with the percentage of ‘A’s, ‘C’s, ‘G’s, ‘T’s, and ‘N’s)
    print("\t%d bases processed (%d%% A, %d%% C, %d%% G, %d%% T, %d%% N)" % (stats["bases"], stats["A"], stats["C"], stats["G"], stats["T"], stats["N"]))
    # total number of bases trimmed (together with the percentage of ‘A’s, ‘C’s, ‘G’s, ‘T’s, and ‘N’s)
    if args.operation == "trim":
        print("\t%d bases trimmed (%d%% A, %d%% C, %d%% G, %d%% T, %d%% N)" % (trimmed_stats["bases"], trimmed_stats["A"], trimmed_stats["C"], trimmed_stats["G"], trimmed_stats["T"], trimmed_stats["N"]))
    # total number of adaptors found (and removed)
    if args.operation == "adaptor-removal":
        print("\t", n_adaptor, "adaptors found")

complement = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}  # complement dictionary
def reverse_complement(sequence):
    rc_seq = ""
    for x in sequence[::-1]:
        rc_seq += complement[x]
    return rc_seq

def trim_right(trim_right):
    right_trimmed_seq = sequence[: -int(trim_right)]
    return right_trimmed_seq

def trim_right_qualities(qualities):
    right_trimmed_qual = qualities[: -int(args.trim_right)]
    return right_trimmed_qual

def trim_left(trim_left):
    left_trimmed_seq = sequence[int(trim_left):]
    return left_trimmed_seq

def trim_left_qualities(qualities):
    left_trimmed_qual = qualities[int(args.trim_left):]
    return left_trimmed_qual

def trim_both(trim_right, trim_left):
    both_trimmed_seq = sequence[int(trim_left):-int(trim_right)]
    return both_trimmed_seq

def trim_both_qualities(qualities):
    both_trimmed_qual = qualities[int(args.trim_left):-int(args.trim_right)]
    return both_trimmed_qual

def remove_adaptor(adaptor):
    global n_adaptor
    adaptor = args.adaptor
    adaptor = adaptor.upper()
    if sequence[0:len(adaptor)] == adaptor:
        n_adaptor += 1
        return sequence[len(adaptor):]
    else:
        return sequence

def check_format(input):
    line = input_file.readline()
    if line[0] == ">":
        input_format = "FASTA"
    elif line[0] == "@":
        input_format = "FASTQ"
    else:
        print("File %s is neither a FASTA/FASTQ file" % input)
        exit(1)
    input_file.seek(0)  # Rewind line
    return input_format  # Return file format


parser = argparse.ArgumentParser(description="FASTX preprocessor")
parser.add_argument("--input", type=str, required=True, help="Input file in FASTA or FASTQ format")
parser.add_argument("--output", type=str, required=True, help="Output file in the same format")
parser.add_argument("--operation", type=str, required=True, choices=["rc", "adaptor-removal", "trim"], help="Choose the operation you would like to perform: reverse-complement(rc), adaptor-removal or trim")
parser.add_argument("--adaptor", type=str, required=False, help="Input an adaptor (e.g.'TATAGA')")
parser.add_argument("--trim-right", type=int, required=False, help="Input length of the trim")
parser.add_argument("--trim-left", type=int, required=False, help="Input length of the trim")
args = parser.parse_args()

if args.operation == "adaptor-removal" and not args.adaptor:
    parser.error("Please input an adaptor (e.g. --adaptor 'TATAGA')")
if args.operation == "adaptor-removal" and args.adaptor.isalpha() == False:
    parser.error("Invalid input. Please input a string (e.g. --adaptor 'TATAGA')")
if args.operation == "trim" and not (args.trim_right or args.trim_left):
    parser.error("Please input whether you want to trim right or left and the number of bases (e.g. --trim-right 10 --trim-left 20)")


# Open input/output files
input_file = open(args.input, "rt")
output_file = open(args.output, "wt")
# Check format
format = check_format(args.input)
# Stats
stats = {"reads": 0, "bases": 0, "A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
trimmed_stats = {"bases": 0, "A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
trimmed_seq = ""
n_adaptor = 0
while True:
    # Read tag or end-of-file
    tag = input_file.readline().rstrip("\n")
    if not tag: break
    # Read sequence
    sequence = input_file.readline().rstrip("\n")
    if args.operation == "adaptor-removal":
        adaptor = args.adaptor
        adaptor = adaptor.upper()
        if len(adaptor) >= len(sequence):
            output_file.write("")
            exit(1)
    if args.operation == "trim":
        if args.trim_right is not None and args.trim_right > len(sequence):
            output_file.write("")
            exit(1)
        elif args.trim_left is not None and args.trim_left > len(sequence):
            output_file.write("")
            exit(1)
        elif (args.trim_right and args.trim_left) is not None and (args.trim_right and args.trim_left) > len(sequence):
            output_file.write("")
            exit(1)
    # FASTQ format
    if format == "FASTQ":
        plusline = input_file.readline().rstrip("\n")
        qualities = input_file.readline().rstrip("\n")
        # Write output file with chosen operation
        if args.operation == "rc":
            output_file.write("%s\n%s\n%s\n%s\n" % (tag, reverse_complement(sequence), plusline, qualities[::-1]))
        elif args.operation == "trim" and args.trim_right:
            output_file.write("%s\n%s\n%s\n%s\n" % (tag, trim_right(args.trim_right), plusline, trim_right_qualities(qualities)))
        elif args.operation == "trim" and args.trim_left:
            output_file.write("%s\n%s\n%s\n%s\n" % (tag, trim_left(args.trim_left), plusline, trim_left_qualities(qualities)))
        elif args.operation == "trim" and args.trim_right and args.trim_left:
            output_file.write("%s\n%s\n%s\n%s\n" % (tag, trim_both(args.trim_right, args.trim_left), plusline, trim_both_qualities(qualities)))
        elif args.operation == "adaptor-removal":
            output_file.write("%s\n%s\n%s\n%s\n" % (tag, remove_adaptor(adaptor), plusline, qualities))
        else:
            print("Operation %s not recognized" % args.operation)
            exit(1)
    # FASTA format
    else:
        # Write output file with chosen operation
        if args.operation == "rc":
            output_file.write("%s\n%s\n" % (tag, reverse_complement(sequence)))
        elif args.operation == "trim" and args.trim_left:
            output_file.write("%s\n%s\n" % (tag, trim_left(args.trim_left)))
        elif args.operation == "trim" and args.trim_right:
            output_file.write("%s\n%s\n" % (tag, trim_right(args.trim_right)))
        elif args.operation == "trim" and args.trim_right and args.trim_left:
            output_file.write("%s\n%s\n" % (tag, trim_both(args.trim_right, args.trim_left)))
        elif args.operation == "adaptor-removal":
            output_file.write("%s\n%s\n" % (tag, remove_adaptor(adaptor)))
        else:
            print("Operation %s not recognized" % args.operation)
            exit(1)
    # Update stats
    stats_update(stats, trimmed_stats, sequence, trimmed_seq)
# Print summary
print_summary(stats, trimmed_stats)
# Close files
input_file.close()
output_file.close()
