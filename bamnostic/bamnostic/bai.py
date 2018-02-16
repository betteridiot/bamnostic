class Bai:
    """ This class defines the bam index file object.
    
    The purpose of this class is the binary parsing of the bam index file (BAI) associated with 
    a given bam file. When queried, the Bai object identifies the bins of data that overlap the 
    requested region and directs which parts of the bam file contain it.
    
    Virtual offsets are processed using the following method:
    Beginning of compressed block = coffset = virtual offset >> 16
    Position within uncompressed block = uoffset = virtual offset ^ (coffset << 16)
    """
    def __init__(self, filename):
        """ Initialization method
        
        Generates an "index" of the index. This gives us the byte positions of each chromosome 
        within the index file. Now, when a user queries over a specific chromosome, it pulls 
        out just the index information for that chromosome--not the whole genome.
        
        Args:
            filename (str): '/path/to/bam_file' that automatically adds the '.bai' suffix
        
        Raises:
            OSError (Exception): if the BAI file is not found or does not exist
        """
        if os.path.isfile(f'{filename}.bai'):
            self.io = open(f'{filename}.bai', 'rb')
        else:
            raise OSError('f{filename}.bai not found. Please change your directory or index your BAM file')
        
        # the 'magic' attribute is just a 4 character string that says the file type in the beginning of the file
        # n_refs is an 'int32_t' type object
        self.magic, self.n_refs = struct.unpack("<4si", self.io.read(struct.calcsize("<4si")))
        
        # the refID_pos attribute is a dictionary of each reference sequence and the byte positions of its corresponding data within the BAI file
        self.refID_pos = {ref[0]:(ref[1], ref[2]) for ref in self.refs(skip = True)}
    
    def offsets(self, n_int):
        """ Pulls out the 16kbp linear indices for a given bin
        
        Args:
            n_int (int): number of intervals listed in the BAI file
            
        Yields:
            ioffset (int): the virtual offset of the first read's position within a given BGZF chunk separated by 16kbp intervals
        """
        for i in range(n_int):
            # each offset is a 'unit65_t' type object
            yield struct.unpack("<Q", self.io.read(struct.calcsize("<Q")))[0]
            
    def chunks(self, n_chunks):
        """ Generator that only pulls out the virtual offsets for the chunk start and stop for all chunks for a given bin
        
        Args:
            n_chunks (int): number of chunks within that given bin
        
        Yields:
            virtual_offsets (tuple): the beginning and end of the a given chunk in virtual offset form
        """
        for c in range(n_chunks):
            chunk_beg, chunk_end = struct.unpack("<2Q", self.io.read(struct.calcsize("<2Q")))
            yield (chunk_beg, chunk_end)   
            
    def bins(self, n_bins, skip = False):
        """ Processes the bin structure of a BAI file
        
        If skip is invoked, then the chunk data is skipped after the number of chunks is determined. Else, the chunk
        data is yielded from the chunks generator
        
        Args:
            n_bins (int): number of bins within the given reference
            skip (bool): whether or not to save the underlying data. This is used to first "index" the index file.
        
        Yields:
            bin_data (dict): bin_id and chunk data key: value pair (if skip is False)
            or
            io.tell() (int): byte position after skipping the chunks (if skip is True)
        """
        for b in range(n_bins):
            bin_id, n_chunk = struct.unpack("<Ii", self.io.read(struct.calcsize("<Ii")))
            
            if skip:
                # skip the chunks data
                self.io.seek(self.io.tell() + (n_chunk * struct.calcsize("<2Q"))) 
            else:
                # delegate to chunks generator for chunk data
                #chunks = yield from self.chunks(n_chunk)
                yield {bin_id: [chunk for chunk in self.chunks(n_chunk)]}
        if skip:
            # send byte position within file after skipping chunks
            yield self.io.tell()
        
    def lookup(self, ref, start, stop):
        """ Main query function for yielding AlignedRead objects from specified region
        
        Args:
            ref (str): which reference/chromosome 
            start (int): left most bp position of region (zero-based)
            stop (int): right most bp position of region (zero-based)
        
        Returns:
            list: offsets of data blocks within BAM file that contain reads that are within specified region
        """
        pass
    
    def refs(self, skip = False, **kwargs):
        """ Main driver function for processing the BAI file
        
        Pulls out the number of bins iteratively and pipelines the pertinent generators to skip or store the underlying data.
        When the skip flag is set to True (default is False), an "index" of the index file is created.
        
        Args:
            skip (bool): Flag for setting whether or not the data is stored
        
        Yields:
            dict: reference id with bin and linear index key: value pair (if skip is False) or (int, int, int) reference id and its byte position
                  start & stop positions within the index (if skip is True)
        """
        for r in range(self.n_refs):
            # n_bin is and 'int32_t' type object 
            n_bin = struct.unpack("<i", self.io.read(struct.calcsize("<i")))[0]
            if skip:
                # get current byte position within file - the n_bin data to capture start of reference block
                ref_bin_start = self.io.tell() - struct.calcsize("<i")
                
                # process the bins while skipping chunks and linear indices
                bins = next(self.bins(n_bin, skip = True))
                n_int = struct.unpack("<i", self.io.read(struct.calcsize("<i")))[0]
                
                # capture the end of the reference block after skipping the linear indices
                ref_bin_end = self.io.seek(self.io.tell() + n_int * struct.calcsize("<Q"))
                yield (r, ref_bin_start, ref_bin_end)
            else:
                # TODO: accept search parameter to isolate reference
                bins = yield from self.bins(n_bin)
                n_int = struct.unpack("<i", self.io.read(struct.calcsize("<i")))[0]
                ioffset = [offset for offset in self.offsets(n_int)]
                yield {r: (bins, ioffset)}