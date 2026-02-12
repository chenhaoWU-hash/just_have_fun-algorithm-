#  SequencesToIndexedSPSS.py -i sequences_file_name -t solidity_threshold -k

import time
import getopt
import sys
import spss
import fmi
import tools_karkkainen_sanders as tks #to sort the suffixes of the sequence
import pickle    #to serialize and deserialize objects for saving FM indexes to files.
import random

## Parse command line arguments

def parse_arguments():
    """
    Analyze command-line arguments with getopt.

    Returns:
        dict: A dictionary containing the parsed arguments.
    
    Raises:
        SystemExit: If required arguments are missing or invalid values are provided.
    """
    # Define short and long options
    short_opts = "hi:k:t:o:"  # h : help (no argument required), i : input, k : kmer_size, t : threshold, o : output 
    long_opts = ["help", "input=", "kmer_size=", "threshold=", "output="]

    # Initialise the parameter dictionary
    args = {
        "input": None,
        "kmer_size": None,
        "threshold": None,
        "output": None
    }

    # Retrieve command line arguments
    argv = sys.argv[1:]

    try:
        opts, _ = getopt.getopt(argv, short_opts, long_opts)
    except getopt.GetoptError as err:
        print(f"Error: {err}")
        print("Use -h or --help for usage information.")
        sys.exit(2)

    for opt, val in opts:
        if opt in ("-h", "--help"):
            print_help()
            sys.exit()
        elif opt in ("-i", "--input"):
            args["input"] = val
        elif opt in ("-k", "--kmer_size"):
            args["kmer_size"] = val
        elif opt in ("-t", "--threshold"):
            args["threshold"] = val
        elif opt in ("-o", "--output"):
            args["output"] = val
        else:
            print(f"Error: Unknown option: {opt}")
            print("Use -h or --help for usage information.")
            sys.exit(2)

    # Check if required parameters are provided
    missing_args = []
    for key in ["input", "kmer_size", "threshold", "output"]:
        if not args[key]:
            missing_args.append(key)
    if missing_args:
        print(f"Error: Missing required arguments: {', '.join(missing_args)}")
        print("Use -h or --help for usage information.")
        sys.exit(2)

    # Input validation
    try:
        args["kmer_size"] = int(args["kmer_size"])
        if args["kmer_size"] <= 0:
            raise ValueError
    except ValueError:
        sys.exit("Error: k-mer size must be a positive integer.")

    try:
        args["threshold"] = int(args["threshold"])
        if args["threshold"] <= 0:
            raise ValueError
    except ValueError:
        sys.exit("Error: Threshold must be a positive integer.")

    return args

def print_help():
    """
    Display the help message for using the script.
    """
    help_message = """
Usage: python SequencesToIndexedSPSS.py -i <input_file> -k <kmer_size> -t <threshold> -o <output_file>

Options:
  -h, --help            Show this help message and exit.
  -i, --input           Input FASTA file containing sequences.
  -k, --kmer_size       Size of k-mers to extract (positive integer).
  -t, --threshold       Solidity threshold for k-mers (positive integer).
  -o, --output          Output file name for the serialized FM-index.

Description:
  Generate SPSS(K) and index it using FM-index.
    """
    print(help_message.strip())

def output(output_file, fm_index):
    """
    Save the FM-index object to a file.

    Args:
        output_file (str): The name of the file to save the FM-index.
        fm_index (object): The FM-index object to serialize and save.
    """
    with open(output_file, "wb") as file:
        pickle.dump(fm_index, file) # Serialize the object and save it to the file

import pickle

def load_fm_index(filename):
    """
    Load the FM index from a file.
    Args:
        filename (str): The name of the file from which to load the FM index.

    Returns:
        object: The deserialized FM-index object.

    Raises:
        FileNotFoundError: If the specified file does not exist.
    """
    try:
        with open(filename, "rb") as file:  # The "rb" mode opens the file in binary read mode
            return pickle.load(file)  # Load the object from the file
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: File {filename} not found.")

def test_fm_index( fm_index, spss_sequence, num_tests=10, k=31):
    """
    Test the validity of the FM-index by comparing the results of the `contains` method
    with Python's `in` operator, using the consensus sequences.

    """
    print("Starting FM-index validation test based on consensus sequences...\n")

    if len(spss_sequence) < k:
        raise ValueError(f"SPSS sequence length ({len(spss_sequence)}) is smaller than k ({k}). Cannot generate valid patterns.")

    success = True
    for i in range(num_tests):
        # Generate a random pattern of length k from the consensus sequence
        start = random.randint(0, len(spss_sequence) - k)
        pattern = spss_sequence[start:start + k]

        fm_result = fm_index.contains(pattern)
       
        # Test with Python `in` on the consensus sequence
        python_result = pattern in spss_sequence

        # Compare results
        if fm_result != python_result:
            print(f"Test failed for pattern '{pattern}' :")
            print(f"FM-index: {fm_result}, Python `in`: {python_result}")
            success = False
        else:
            print(f"Test passed for pattern '{pattern}'.")

    if success:
        print("\nAll tests passed successfully.")
    else:
        print("\nDiscrepancies were detected between the FM-index and Python `in`.")
    return success



def main():
    """
    Main function to orchestrate the script operations.

    Steps:
        1. Parse command-line arguments.
        2. Extract canonical k-mers from input sequences.
        3. Filter k-mers based on solidity threshold.
        4. Extend k-mers to generate SPSS(K).
        5. Concatenate SPSS(K) and build an FM-index.
        6. Serialize and save the FM-index to the specified output file.
    """
    args = parse_arguments()
    print("Program is running...")
    
    # Retrieve parameters
    input_file = args['input']
    kmer_size = args['kmer_size']
    threshold = args['threshold']
    output_file = args['output']

    t0 = time.time()
    
    #  Build SPSS
    canonical_kmers = spss.listing_kmer(input_file, kmer_size)
    solide_canonical_kmers = spss.solide(canonical_kmers, threshold)
    t1 = time.time()
    
    spss_list = spss.extend(solide_canonical_kmers)
    spss_sequence = spss.concat_spss(spss_list)  # Consensus SPSS
    t2 = time.time()

    #  Build the FM-index
    my_fm_index = fmi.FmIndex(spss_sequence)
    bwt = my_fm_index.set_bwt()
    n, rank = my_fm_index.set_n_and_rank()
    t3 = time.time()

    # Display timing statistics
    print(f"OUT TIME_SELECTING_KMERS={t1-t0}")
    print(f"OUT |SPSS(K)|={len(''.join(spss_list))}")
    print(f"OUT #SPSS(K)={len(spss_list)}")
    print(f"OUT TIME_SPSS_CONSTRUCTION={t2-t1}")
    print(f"OUT TIME BUILD FMI={t3-t2}")

    # Save the FM-index
    output(output_file, my_fm_index)

    # Validate the FM-index
    print("\nValidating FM-index...")
    try:
        # Load the raw sequence
        with open(input_file, 'r') as f:
            original_sequence = "".join(line.strip() for line in f if not line.startswith(">"))
        
        # Test the FM-index with consensus sequences
        test_result = test_fm_index(my_fm_index, spss_sequence)
        
        if test_result:
            print("Validation successful: FM-index matches Python `in` operator for all tests.")
        else:
            print("Validation failed: Discrepancies found between FM-index and Python `in` operator.")
    except Exception as e:
        print(f"An error occurred during FM-index validation: {e}")


if __name__ == "__main__":
    main()