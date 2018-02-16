class BAMheader:
    def __init__(self, reader):
        self.magic, self.header_length = struct.unpack('<4si', reader.read(8))
        
        # If SAM header is present, it is in plain text. Process it and save it as rows
        self.SAMheader = struct.unpack(f'<{self.header_length}s', reader.read(self.header_length))[0].decode()
        self.SAMheader = [i.split("\t") for i in self.SAMheader.split("\n")]
        self.SAMheader_end = reader.tell()
        
        # Each reference is listed with the @SQ tag. We need the number of refs to process the data
        self.n_refs = struct.unpack('<i', reader.read(4))[0]
        
        # create a dictionary of all the references and their lengths
        self.refs = {}
        for r in range(self.n_refs):
            name_len = struct.unpack("<i", reader.read(4))[0]
            ref_name = struct.unpack(f'{name_len-1}s', reader.read(name_len)[:-1])[0]
            ref_len = struct.unpack("<i", reader.read(4))[0]
            self.refs.update({r: (ref_name.decode(), ref_len)})
        self.BAMheader_end = reader.tell()
        print(f'End position of BAM header: {self.BAMheader_end}')
    
class AlignedRead:
    """Main class for handling reads within the BAM"""
    
    def __init__(self, reader):
        """Instantiating the read parser just needs access to the BGZF reader object
        
        Args:
            reader (BgzfReader): parser for processing BGZF files
            
        Returns:
            AlignedRead
        
        """
        self.block_size = struct.unpack('<i', reader.read(4))[0]
        self.byte_stream = bytearray(reader.read(self.block_size))
        self.refID, self.pos = struct.unpack('<2i', self.range_popper(2*4))
        
        # get refID chromosome names
        self.refID = header.refs[self.refID][0]
        self.bin_mq_nl = struct.unpack('<I', self.range_popper(4))[0]
        self.bin = self.bin_mq_nl >> 16
        self.mapq = (self.bin_mq_nl & 0xFF00) >> 8
        self.l_read_name = self.bin_mq_nl & 0xFF
        
        # Alternative to masking
        # self.mapq = (self.bin_mq_nl ^ self.bin << 16) >> 8
        # self.l_read_name = (self.bin_mq_nl ^ self.bin <<16) ^ (self.mapq << 8)
        
        self.flag_nc = struct.unpack('<I', self.range_popper(4))[0]
        self.flag = self.flag_nc >> 16
        self.n_cigar_op = self.flag_nc & 0xFFFF
        self.l_seq, self.next_refID, self.next_pos, self.tlen = struct.unpack('<4i', self.range_popper(4*4))
        self.read_name = struct.unpack(f'<{self.l_read_name}s', self.range_popper(self.l_read_name))[0].decode()[:-1]
        
        # Both the CIGAR and seq strings require determining size, decoding, and mapping to key
        cigar_key = "MIDNSHP=X"
        self.cigar = struct.unpack(f'<{self.n_cigar_op}I', self.range_popper(4 * self.n_cigar_op))
        self.cigar_string = "".join([f'{self.cigar[i]>> 4}{cigar_key[self.cigar[i] & 0xF]}' for i in range(len(self.cigar))])
        
        self.byte_seq = struct.unpack(f'<{(self.l_seq + 1)//2}B', self.range_popper(1 * ((self.l_seq + 1)//2)))
        
        seq_key = '=ACMGRSVTWYHKDBN'
        self.seq = "".join([f'{seq_key[self.byte_seq[s] >> 4]}{seq_key[self.byte_seq[s] & 0x0F]}' for s in range(len(self.byte_seq))])[:self.l_seq]

        
        self.qual = struct.unpack(f'<{self.l_seq}s', self.range_popper(self.l_seq))[0]
        self.tags = {}
        while len(self.byte_stream) > 0:
            self.tags.update(self.tagger())
        
    def range_popper(self, rng):
        """Simple pop method that accepts a range instead of a single value. Modifies the original bytearray by removing items
        
        Args:
            rng (int): desired number of bytes from the beginning to be decoded
            
        Returns:
            popped (bytearray): removed rng-length items from byte_stream 
        
        """
        popped = self.byte_stream[:rng]
        del self.byte_stream[:rng]
        return popped
    
    def tagger(self):
        """Processes the variety of tags that accompany sequencing reads
        
        Returns:
            dictionary of the tag, value types, and values
        """
        types = {"A":'<c',"i":'<l',"f":'<f',"Z":'<s',"H":'<s',"c":'<b',"C":'<B',"s":'<h',"S":'<H',"i":'<i',"I":'<I'}
        tag = struct.unpack('<2s', self.range_popper(2))[0].decode()
        val_type = struct.unpack('<s', self.range_popper(1))[0].decode()
        if val_type == "B":
            arr_type = struct.unpack('<s', self.range_popper(1))[0].decode()
            arr_size = struct.unpack('<i', self.range_popper(4))[0]
            arr = struct.unpack(f'<{arr_size}{types[arr_type]}', self.range_popper(arr_size * struct.calcsize(types[arr_type])))
            return {tag: (val_type, arr)}
        elif val_type == "Z" or val_type == "H":
            val = struct.unpack('<s', self.range_popper(1))[0]
            while val[len(val)-1:len(val)] != b'\0':
                val += struct.unpack('<s', self.range_popper(1))[0]
            if val_type == "Z":
                return {tag: (val_type, val.decode(encoding="ascii")[:-1])}
            else:
                return {tag: (val_type, val.hex().upper())}
        else:
            val_size = struct.calcsize(types[val_type])
            val = struct.unpack(types[val_type], self.range_popper(val_size))[0]
            if val_type == "A":
                val = val.decode(encoding='ascii')
            return {tag: (val_type, val)}
        