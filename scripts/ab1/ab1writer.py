#!/usr/bin/env python3

"""
Adapted from Brett Milash's script (https://github.com/bmilash)
"""

### ---------------------------------------- ###

def FreqToBase( base_frequencies ):
    """Given a list of base frequency dictionaries, return a string with
    the most abundant base at each position."""
    return ''.join([ max(x.items(),key=lambda y:y[1])[0] for x in base_frequencies ])

def CalculateNucleotideFrequencies(sequence):
    """Returns dictionary of nucleotides and their relative frequencies
    in a sequence. For now ignoring insertions and deletions."""
    return { i:sequence.upper().count(i)/len(sequence) for i in 'ACGT' }

def PileupToBaseFrequencyList(pileup_filename):
    """Given the name of a samtools pileup file, returns list of
    nucleotide frequency dictionaries, one dictionary for each
    position of the pileup."""
    base_frequency_list=[]
    with open(pileup_filename) as ifs:
        for line in ifs:
            seqdata=line.split('\t')[4]
            base_frequency_list.append(CalculateNucleotideFrequencies(seqdata))
    return base_frequency_list

def StringToBaseFrequencyList(seqdata):
    """Given a DNA sequence as a string, returns a list of
    nucleotide frequency dictionaries, one dictionary for each
    position of the sequence."""
    base_frequency_list=[]
    for nucleotide in seqdata:
        base_frequency_list.append(CalculateNucleotideFrequencies(seqdata))
    return base_frequency_list

def GenerateFakeQualities(sequence):
    """Creates a fake quality score sequence - just a string of H's as long as the sequence."""
    return 'H' * len(sequence)

def GenerateFakeIntensities(base_frequencies,nucleotide,samples=4):
    """
    Intend to use 4 samples per base looking up intensities from a 
    dictinary.
    """
    intensities={
        15:[ 0, 177, 678, 1415, 2262, 3072, 3705, 4051, 4051, 3705, 3072, 2262, 1415, 678, 177 ],
        4:[ 0, 2048, 4096, 2048 ],
    }
    i_values=np.array(intensities[samples])
    # Create 2D array for intensity data. One row per nucleotide in the
    # sequence, and one column per intensity value. Using 15 intensities
    # per nucleotide to creates nice curved peaks in the trace, but this
    # limits the length of the sequences we can store, so backing off to
    # 4 intensities per nucleotide.
    intensities = np.zeros((len(base_frequencies),samples),dtype=np.int16)
    for (i,freq) in enumerate(base_frequencies):
        intensities[i,:]=i_values*freq[nucleotide]
    return intensities.reshape((len(base_frequencies)*samples,))

class AbiHeader:
    """AbiHeader represents the header structure for an ABIF file."""
    header_format=">4sH4sI2H3I"
    #header_size=struct.calcsize(header_format)
    # Fixing header size at 128 bytes.
    # Header size is 28 bytes, but the area of the file is padded
    # out to 128 bytes.
    header_size=28

    @staticmethod
    def read(ifs):
        """Reads from file object opened in "rb" mode and returns an AbiHeader object."""
        header_values=struct.unpack(AbiHeader.header_format, ifs.read(struct.calcsize(AbiHeader.header_format)))
        return AbiHeader( *header_values )

    #(101, b'tdir', 1, 1023, 28, 130, 4032, 205192)
    def __init__(self, abif_tag=b'ABIF', x1=101, x2=b'tdir', x3=1, x4=1023, header_size=header_size, num_entries=0, x7=4032, directory_location=header_size):
        self.abif_tag=abif_tag
        self.x1=x1
        self.x2=x2
        self.x3=x3
        self.x4=x4
        self.header_size=header_size
        self.num_entries=num_entries
        self.x7=x7
        self.directory_location=directory_location

    def write(self,ofs):
        """Writes the ABIF file header to the output file stream, which must
        be opened in 'wb' mode."""
        ofs.write(struct.pack( self.header_format,
            self.abif_tag,
            self.x1,
            self.x2,
            self.x3,
            self.x4,
            self.header_size,
            self.num_entries,
            self.x7,
            self.directory_location))
        # Pad the file to the full size of the header.
        #padsize=self.header_size - struct.calcsize(self.header_format)
        padsize=98
        pad_format=f'>{padsize}s'
        ofs.write( struct.pack( pad_format, bytes([0] * padsize )))
    
    def __repr__(self):
        return f"abif_tag={self.abif_tag},x1={self.x1},x2={self.x2},x3={self.x3},x4={self.x4},header_size={self.header_size},num_entries={self.num_entries},x7={self.x7},directory_location={self.directory_location}"

class DataBlock:
    """
    DataBlock is an abstract class that helps identify data block objects
    as such (as opposed to directory entry objects.
    """
    pass

class DirectoryEntry:
    """DirectoryEntry represents one Directory Entry at the end of
    an ABIF file. The structure can either hold data in-line or
    point to the location of the data in the ABIF file."""
    dir_format = ">4sI2H4I"
    table_header='\t'.join(["name","number","type","el_size","num_el","data","offset","handle"])

    @staticmethod
    def read(ifs):
        """Reads from file object opened in "rb" mode and returns a DirectoryEntry object."""
        dir_entry=struct.unpack(DirectoryEntry.dir_format, ifs.read(struct.calcsize(DirectoryEntry.dir_format)))
        return DirectoryEntry( *dir_entry )

    def __init__(self,name,number,elementtype,elementsize,numelements,datasize,dataoffset,datahandle):
        self.name=name
        self.number=number
        self.elementtype=elementtype
        self.elementsize=elementsize
        self.numelements=numelements
        self.datasize=datasize
        self.dataoffset=dataoffset
        self.datahandle=datahandle
    
    def __lt__(self,other):
        return self.sort_by_name(other)

    def sort_by_location(self,other):
        return self.dataoffset < other.dataoffset

    def sort_by_name(self,other):
        if self.name == other.name:
            return self.number < other.number
        return self.name < other.name

    def __repr__(self):
        return f"{self.name}\t{self.number}\t{self.elementtype}\t{self.elementsize}\t{self.numelements}\t{self.datasize}\t{self.dataoffset}\t{self.datahandle}"

    def write(self,ofs):
        """
        Writes directory entry to file object opened in write binary mode.
        """
        ofs.write(struct.pack( self.dir_format,
            self.name,
            self.number,
            self.elementtype,
            self.elementsize,
            self.numelements,
            self.datasize,
            self.dataoffset,
            self.datahandle ))

class DATA(DataBlock):
    """
    DATA records store the intensity data for each channel.
    name    number    type    el_size    num_el        data        offset    handle
    b'DATA'    9    4    2    4*seqlen    8*seqlen    180536    0
    b'DATA'    10    4    2    4*seqlen    8*seqlen    180536    0
    b'DATA'    11    4    2    4*seqlen    8*seqlen    180536    0
    b'DATA'    12    4    2    4*seqlen    8*seqlen    180536    0
    """
    def __init__(self,number,data):
        self.data=data
        self.dir_entry=DirectoryEntry(name=b'DATA',
            number=number,
            elementtype=4,
            elementsize=2,
            numelements=len(self.data),
            datasize=len(self.data)*2,
            dataoffset=None,
            datahandle=0)
        
    def set_offset(self,offset_value):
        self.dir_entry.dataoffset = offset_value

    def write(self,ofs):
        """
        Writes DATA record to file, and records the offset where
        the entry was written.
        """
        self.dir_entry.offset=ofs.tell()
        packing_format=f">{len(self.data)}H"
        ofs.write(struct.pack(packing_format,*self.data))

class PCON(DataBlock):
    """
    PCON records store the per-base quality score.
    name    number    type    el_size    num_el    data    offset    handle
    b'PCON'    1    2    1    795    795    180536    0
    b'PCON'    2    2    1    795    795    180536    0
    """
    def __init__(self,qualities,number):
        self.data=qualities
        self.dir_entry=DirectoryEntry(name=b'PCON',
            number=number,
            elementtype=2,
            elementsize=1,
            numelements=len(qualities),
            datasize=len(qualities),
            dataoffset=None,
            datahandle=0)
        
    def set_offset(self,offset_value):
        self.dir_entry.dataoffset = offset_value

    def write(self,ofs):
        """
        Writes PCON record to file, and records the offset where
        the entry was written.
        """
        self.dir_entry.offset=ofs.tell()
        ofs.write(self.data.encode('utf-8'))

class PLOC(DataBlock):
    """
    PLOC stores the locations of the peaks in the trace. Its just
    an array of shorts.
    name    number    type    el_size    num_el    data    offset    handle
    b'PLOC'    1    4    2    5417    10834    5545    0
    """
    def _Data1( self, sequence ):
        # Generate the data for the PLOC data block. This needs to
        # blocks of 64 values, which are multiples of 1024 + an offset.
        # The offset is 768 + the block number (which starts at 0).
        # Doing this as a 2D numpy array, 64 columns, and the number
        # of rows is len(sequence)//64 + 1.
        numblocks=len(sequence)//64+1
        numcols=64
        #data=np.zeros((numblocks,numcols),dtype=np.uint16)
        data=np.zeros((numblocks,numcols)).astype(np.uint16)

        #data=np.fromfunction(lambda x,y:y*1024+x,(3,64),dtype=np.uint16) + 768
        data=np.fromfunction(lambda x,y:y*1024+x,(3,64)).astype(np.uint16) + 768
        self.data=data.flatten()

    def _Data2(self,sequence):
        return [ i*4+3 for i in range(len(sequence)) ]

    def _Data3(self,sequence,samples=4):
        """Data3 returns peak locations with 4 samples (default)
        per nucleotide."""
        #data = np.array([ i*samples+(samples//2) for i in range(len(sequence)) ],dtype=np.uint16)
        data = np.array([ i*samples+(samples//2) for i in range(len(sequence)) ]).astype(np.uint16)
        assert(len(data) == len(sequence))
        return data

    def __init__(self,sequence,number=1):
        self.data=self._Data3(sequence)
        self.dir_entry=DirectoryEntry(name=b'PLOC',
            number=number,
            elementtype=4,
            elementsize=2,
            numelements=len(self.data),
            datasize=len(self.data)*2,
            dataoffset=None,
            datahandle=0)
        
    def set_offset(self,offset_value):
        self.dir_entry.dataoffset = offset_value

    def write(self,ofs):
        """
        Writes PLOC record to file, and records the offset where
        the entry was written.
        Note - can only write about 8,000 peak locations at a
        time using struct.pack, so need to break up data into
        chunks <= 8k in length.
        """
        self.dir_entry.offset=ofs.tell()
        chunksize=100
        start=0
        end=min(chunksize,self.dir_entry.numelements)
        while start < self.dir_entry.numelements:
            
            format_string=f">{end-start}H"
            #print(f"start {start} end {end} format_string {format_string}")
            ofs.write(struct.pack(format_string,*list(self.data[start:end])))
            start=end
            end=min(end+chunksize,self.dir_entry.numelements)


class PBAS(DataBlock):
    """
    PBAS records store the base-called sequence data.
    name    number    type    el_size    num_el    data    offset    handle
    b'PBAS'    1    2    1    795    795    180536    0
    b'PBAS'    2    2    1    795    795    180536    0
    """
    def __init__(self,sequence,number):
        self.data=FreqToBase(sequence)
        self.dir_entry=DirectoryEntry(name=b'PBAS',
            number=number,
            elementtype=2,
            elementsize=1,
            numelements=len(sequence),
            datasize=len(sequence),
            dataoffset=None,
            datahandle=0)
        
    def set_offset(self,offset_value):
        self.dir_entry.dataoffset = offset_value

    def write(self,ofs):
        """
        Writes PBAS2 record to file, and records the offset where
        the entry was written.
        """
        self.dir_entry.offset=ofs.tell()
        ofs.write(self.data.encode('utf-8'))

class FWO(DirectoryEntry):
    """
    FWO is a kind of DirectoryEntry record that stores the base order,
    ie which base is Channel 1, which base is Channel 2, and so on.
    The FWO directory entry holds the four bases, in order, in a 32-bit
    int, and this information is stored within the FWO entry in the
    directory (in the dataoffset field) rather than storing an offset to 
    the information in a data block the file.
    """
    @staticmethod
    def Encode(bases):
        """
        Translates character string of 4 nucleotides (e.g. 'ACGT') 
        into 32-bit int holding the ascii values of those characters.
        """
        encoded_value=0
        for base in bases:
            encoded_value = (encoded_value << 8)+ord(base)
        return encoded_value

    @staticmethod
    def Decode(encoded_value):
        """
        Translates 32-bit int back to a 4-character string with
        the 4 nucleotides in the order in which they are stored
        in the ABIF file.
        """
        bases=[]
        for i in range(4):
            bases.append( chr(encoded_value & 0xff))
            encoded_value = ( encoded_value >> 8 )
        bases.reverse()
        return ''.join(bases)
    
    def __init__(self,bases):
        super().__init__(name=b'FWO_',
            number=1,
            elementtype=2,
            elementsize=1,
            numelements=4,
            datasize=4,
            dataoffset=FWO.Encode(bases),
            datahandle=0)
        self.channel_mapping = { bases[0]:9, bases[1]:10, bases[2]:11, bases[3]:12 }

    def __getitem__(self,key):
        return self.channel_mapping[key]
    
class Trace:
    """
    Trace represents a single ABI file trace which stores sequence and quality data. The object's
    constructor takes a single required argument: the DNA sequence data. Additional optional arguments
    control resolution, base order, etc. The single public method, write() writes the Trace object to
    a .ab1 file.
    One could use the Trace object like this:
    Trace(seqdata).write(output_file_name)
    """
    def __init__(self,seqdata,baseorder='GATC',resolution=15,pileup=True):
        # Need to initialize the following data blocks:
        # name    number
        # b'DATA'    9
        # b'DATA'    10
        # b'DATA'    11
        # b'DATA'    12
        # b'FWO_'    1
        # b'PBAS'    1
        # b'PBAS'    2
        # b'PCON'    1
        # b'PCON'    2
        # b'PLOC'    1

        # Transform the sequence data into a list of dictionaries
        # with the relative frequency of each nucleotide at each
        # position.
        if pileup:
            self.base_frequencies=PileupToBaseFrequencyList(seqdata)
        else:
            self.base_frequencies=StringToBaseFrequencyList(seqdata)
        
        # Create the file header.
        self.abi_header = AbiHeader()

        # Create the base order object.
        self.fwo=FWO(baseorder)
        self.data_blocks= [ PBAS(self.base_frequencies,1), 
            PLOC(self.base_frequencies),
            PBAS(self.base_frequencies,2), 
            ]
        # Store the trace length which is the length of the PBAS data.
        self.trace_length = len(self.data_blocks[0].data)
        # Store the consensus sequence.
        self.consensus_sequence=self.data_blocks[0].data
        # Create the DATA and PCON blocks.
        for base in baseorder:
            self.data_blocks.append(DATA(self.fwo[base],GenerateFakeIntensities(self.base_frequencies,base)))
        self.data_blocks.append(PCON(GenerateFakeQualities(self.base_frequencies),1))
        self.data_blocks.append(PCON(GenerateFakeQualities(self.base_frequencies),2))

        # Initialize the directory with the directory entries that 
        # store their data in-line (rather than point to a data block).
        self.directory = [ self.fwo, ]

        # Calculate position of each data block in file starting 
        # location of first data block - first data written at 
        # offset 128 due to unused bytes following header.
        position = 128

        # For each data block ...
        for block in self.data_blocks:
            # ... record block's location...
            block.set_offset(position)
            # ... and calculate position of subsequent block.
            position+=block.dir_entry.datasize
            # Store data block's directory entry in directory.
            self.directory.append(block.dir_entry)
        # Record the location of the directory (which follows all the data blocks).
        self.abi_header.directory_location=position
        # Record the number of directory entries.
        self.abi_header.num_entries = len(self.directory)
    
    def write(self,filename):
        """
        Writes trace to a .ab1 file.
        """
        with open(filename,'wb') as ofs:
            # Write the header.
            self.abi_header.write(ofs)

            # Write each block.
            for block in self.data_blocks:
                block.write(ofs)

            # Sort and write the directory.
            self.directory.sort()
            for dir_entry in self.directory:
                dir_entry.write(ofs)

### ------------------MAIN------------------ ###

import struct
import string
import itertools
import numpy as np

from sys import argv

# Parse args

input_file = argv[argv.index('--pileup') + 1]

output_file = argv[argv.index('--output') + 1]

# Generate ab1 trace

synthetic_trace = Trace(input_file, pileup=True)

synthetic_trace.write(output_file)
