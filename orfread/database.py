from sqlalchemy import Integer, String, Column, Float, ForeignKey
from sqlalchemy.orm import declarative_base, relationship
import datetime

Base=declarative_base()

#RAW file database table
#this is processed outside by the user
class raw_file(Base):
    __tablename__='raw_files'
    #RAW file database schema
    id=Column(Integer, primary_key=True)
    file_name=Column(String, nullable=False) #file name header
    header=Column(String, nullable=False) #sequencing file header
    orf_proteins=Column(Integer)
    date=Column(String, default=datetime.datetime.now()) #date time of entry

    #relationships
    orfs=relationship('ORF', back_populates='raw_file', #one to many relationship, as one file can have many ORF proteins
                      cascade='all, delete-orphan') #deletes all related ORF proteins, if the raw_file is deleted too

    #represent method to print data
    def __repr__(self):
        return(f'file_name:{self.file_name}, header:{self.header}, date:{self.date}')


#ORF database table
class ORF(Base):
    __tablename__='ORFs'
    #ORF database schema
    ORF_id=Column(Integer, primary_key=True)
    ORF_name=Column(String) #ORF_name eg: 'ORF PROTEIN 1'
    ORF_header=Column(String, ForeignKey('raw_files.id')) #ORF Protein 1 etc... with a relation to raw_files.id (foreign key)
    ORF_strand=Column(String) #ORF strand, forward (5'>'3' Upstream=Forward strand) or reverse (3'<5' Downstream=Reverse strand)
    ORF_frameshift=Column(Integer) #ORF Frame shift of upstream or downstream shifts of 1, 2 or 3
    #ORF_frame_shift=Column(String)#ORF Frame shift eg: 5'>3' Frame 1, 5'>3' Frame 2, 3'<5' Frame 1...etc MIGHT DELETE THIS IF IM GOING WITH ORF_strand and ORF_frameshift

    ORF_seq=Column(String) #ORF sequence, eg: MASSLAJWAL-
    ORF_index=Column(Integer) #ORF start codon 'M' index position on TRANSLATED file, if i delete translated files, i will have to delete this as well
    ORF_len=Column(Integer) #ORF sequence length

    NTD_seq=Column(String) #RAW (true) nucleotide sequence of ORF seq (taken from raw file, meaning its NOT simply translated from ORF protein, so it can contain RAW mutations of RAW seq)
    NTD_index=Column(Integer) #NTD sequence start codon index position in the raw file
    NTD_len=Column(Integer) #NTD sequence length

    #relationships
    raw_file=relationship("raw_file", back_populates="orfs") #many to one relationship, as many ORF proteins belong/come from the same raw sequenced file

    #represent method to print data for testing
    def __repr__(self):
        return(f'ORF_id:{self.ORF_id}, ORF_header:{self.ORF_header}, ORF_strand:{self.ORF_strand}, ORF_frameshift:{self.ORF_frameshift}, ORF_index:{self.ORF_index}, ORF_len:{self.ORF_len}, NTD_index:{self.NTD_index}')


if __name__=='__main__':

    print('hi')