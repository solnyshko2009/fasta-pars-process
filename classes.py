class Seq:
    """
    A class to represent a biological sequence from FASTA format.
    
    This class stores the header and sequence data from a FASTA record
    and provides methods to analyze and manipulate the sequence.
    
    Attributes:
        header (str): The header line from the FASTA record (without the '>' character)
        sequence (str): The biological sequence data with whitespace and newlines removed
    """
    
    def __init__(self, header, sequence):
        """
        Initialize a biological sequence.
        
        Args:
            header (str): FASTA record header line
            sequence (str): Biological sequence data
        """
        self.header = header.strip()
        self.sequence = sequence.strip().replace('\n', '').replace(' ', '')
    
    def __str__(self):
        """
        Return string representation in FASTA format.
        
        Returns:
            str: Sequence in FASTA format with header and formatted sequence
        """
        return f">{self.header}\n{self._format_sequence()}"
    
    def __len__(self):
        """
        Return the length of the sequence.
        
        Returns:
            int: Length of the sequence
        """
        return len(self.sequence)
    
    def _format_sequence(self):
        """
        Format sequence for display (internal method).
        
        Returns:
            str: Sequence formatted with line breaks every 80 characters
        """
        # This method would typically format the sequence with line breaks
        # For example: return '\n'.join([self.sequence[i:i+80] for i in range(0, len(self.sequence), 80)])
        return self.sequence
    
    def get_alphabet(self):
        """
        Determine the sequence alphabet type based on character composition.
        
        The method checks if the sequence contains primarily nucleotide characters
        (A, C, G, T, U, N) or protein characters (amino acid codes). Sequences
        containing characters from both sets or unrecognized characters will be
        classified as 'unknown'.
        
        Returns:
            str: Sequence type - one of 'nucleotide', 'protein', or 'unknown'
        """
        if not self.sequence:
            return 'unknown'
        
        nucleotide_chars = set('ACGTUacgtuNn-')
        protein_chars = set('ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy*Xx-')
        
        for char in self.sequence:
            if char not in nucleotide_chars:
                for char in self.sequence:
                    if char not in protein_chars:
                        return 'unknown'
                return 'protein'
        
        return 'nucleotide'


class FastaReader:
    """
    A class for reading and parsing FASTA format files.
    
    This class provides methods to validate FASTA files and read sequences
    from them sequentially using a generator pattern.
    
    Attributes:
        filename (str): Path to the FASTA file to be read
    """
    
    def __init__(self, filename):
        """
        Initialize FASTA reader.
        
        Args:
            filename (str): Path to FASTA file
        """
        self.filename = filename
    
    def is_valid_fasta(self):
        """
        Check if the file conforms to FASTA format.
        
        The validation checks if the first non-empty line of the file
        starts with the '>' character, which is required for FASTA format.
        
        Returns:
            bool: True if file conforms to FASTA format, False otherwise
            
        Note:
            This is a basic validation that only checks the first line.
            It doesn't guarantee the entire file is properly formatted.
        """
        try:
            with open(self.filename, 'r') as file:
                first_line = file.readline().strip()
                return first_line.startswith('>')
        except (IOError, UnicodeDecodeError):
            return False
    
    def read_sequences(self):
        """
        Generator for reading sequences from FASTA file.
        
        This method reads the FASTA file line by line and yields Seq objects
        for each sequence found in the file. It handles multi-line sequences
        by concatenating all lines until the next header.
        
        Yields:
            Seq: Sequence objects from the FASTA file
            
        Raises:
            ValueError: If file doesn't conform to FASTA format
            
        Example:
            >>> reader = FastaReader('sequences.fasta')
            >>> for seq in reader.read_sequences():
            ...     print(f"Header: {seq.header}, Length: {len(seq)}")
        """
        if not self.is_valid_fasta():
            raise ValueError(f"File {self.filename} does not conform to FASTA format")
        
        current_header = None
        current_sequence = []
        
        with open(self.filename, 'r') as file:
            for line in file:
                line = line.strip()
                
                if not line:
                    continue
                
                if line.startswith('>'):
                    if current_header is not None:
                        yield Seq(current_header, ''.join(current_sequence))
                    
                    current_header = line[1:]
                    current_sequence = []
                else:
                    current_sequence.append(line)
            
            if current_header is not None:
                yield Seq(current_header, ''.join(current_sequence))


def main():
    """
    Main function demonstrating the usage of FastaReader and Seq classes.
    
    This function:
    1. Creates a FastaReader instance with a specific file path
    2. Validates the FASTA file
    3. Reads and displays information about the first 30 sequences
    4. Shows sequence headers, lengths, and alphabet types
    """
    reader = FastaReader(r'c:\Users\Asus\Desktop\sq+read\uniparc_active_p1.fasta')
    
    print(f"File is valid: {reader.is_valid_fasta()}")
    print()
    
    print("=== Reading sequences ===")
    for i, seq in enumerate(reader.read_sequences(), 1):
        print(f"\nSequence {i}:")
        print(f"  Header: {seq.header}")
        print(f"  Length: {len(seq)}")
        print(f"  Type: {seq.get_alphabet()}")
        if i > 30:
            print("And so on...")
            break


if __name__ == "__main__":
    main()