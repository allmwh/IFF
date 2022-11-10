import re
import git
import numpy as np
from Bio import SeqIO

from function.utilities import find_human_sequence


class SeqFilter():
    '''
    some od_ident (derived from PONDR) postprocessing functions 

    symbols used in od_ident：
        1: disorder identified by PONDR
        0: order identified by PONDR

        x: disorder identified by PONDR, but the length is shorter than condition by length_filter_by_od_ident() 
        z:    order identified by PONDR, but the length is shorter than condition by length_filter_by_od_ident() 
    '''
    def __init__(self):
        pass

    def length_filter_by_od_ident(self, od_ident, disorder_filter_length, order_filter_length):
        '''
        filter order/disorder seqeunce which is shorter than specified condition

        od_ident: str
        disorder_filter_length: int, continuous disorder sequence which is shorter than is replaced with "x" 
        order_filter_length: int, continuous order sequence which is shorter than is replaced with "z"

        return: od_ident, str
        '''
        #filter_length
        disorder_check = re.finditer("1+", od_ident)
        for i in disorder_check:
            if i:
                start = i.start()
                end = i.end()
                if end - start <= disorder_filter_length:
                    od_ident = od_ident[:start] + "x" * (end - start) + od_ident[end:]
        
        order_check = re.finditer("0+", od_ident)
        for i in order_check:
            if i:
                start = i.start()
                end = i.end()
                if end - start <= order_filter_length:
                    od_ident = od_ident[:start] + "z" * (end - start) + od_ident[end:]

        return od_ident

    def get_seq_from_od_ident(self, od_ident,sequence,od):
        '''
        get order/disorder masked protein sequence from od_ident 
        the lengths of od_ident and sequence must be same
        
        od_ident: str
        sequence: str, protein seqeunce 
        od: str, ['order', 'disorder'] make order or disorder sequence 
        '''
        new_seq = ''
        for index, element in enumerate(od_ident):
            if od == 'order': # make order seq
                if element == '0':
                    new_seq = new_seq + sequence[index]
                elif element == 'z':
                    new_seq = new_seq + 'z'
                elif element == '1' or element == 'x':
                    new_seq = new_seq + '*'
            
            elif od == 'disorder':# make disorder seq
                if element == '1':
                    new_seq = new_seq + sequence[index]
                elif element == 'x':
                    new_seq = new_seq + 'x'
                elif element == '0' or element == 'z':
                    new_seq = new_seq + '*'
                    
        return new_seq

    def get_od_index(self, od_ident):
        '''
        get order/disorder index 

        od_ident: str

        return: order/disorder index   
        '''
        order_region = []
        disorder_region = []

        disorder_check = re.finditer("1+", od_ident)
        for i in disorder_check:
            if i:
                start = i.start()
                end = i.end()

                disorder_region.append({"start": start, "end": end})

        order_check = re.finditer("0+", od_ident)
        for i in order_check:
            if i:
                start = i.start()
                end = i.end()

                order_region.append({"start": start, "end": end})

        return {"order_region": order_region, 
                "disorder_region": disorder_region}

    def od_add_alignment(self, path, od_ident):
        '''
        add the gap to the od_ident from human protein sequence 
        to ensure that the length between od_ident and homologs sequences from alied fasta file are the same,
        
        there are many conditions are considered:
            
            nen1: 1---1, both sides of the gaps are 1 (disorder idetified by PONDR)
            nen0: 0---0, both sides of the gaps are 0 (order idetified by PONDR)
            nenx: x---x, both sides of the gaps are x (order idetified by PONDR, but the length does not satisfy the filter condition from length_filter_by_od_ident())
            nenz: z---z, both sides of the gaps are z (disorder idetified by PONDR, but the length does not satisfy the filter condition from length_filter_by_od_ident())
            een: ---(0, 1, x, z), sequence starting from the gaps, filling 0 (order) ,1 (disorder), x (disorder filtered) or z (order filtered) respectively 
            nee: (0, 1, x, z)---, sequence ending from the gaps, filling 0 (order) ,1(disorder), x (disorder filtered) or z (order filtered) respectively 
            
            complicated condition:
            nen01: (0, 1)---(0, 1), different identifier (0 or 1) of both sides of the gaps, we filled 0 (order)
            nenxz: (x, z)---(x, z), different identifier (x or z) of both sides of the gaps, we filled z (order filtered)
            
            nen1z: (1, z)---(z, 1), both sides of the gaps are disorder, but one side is z, we filled z (disorder filtered)
            nen0x: (0, x)---(x, 0), both sides of the gaps are order, but one side is x, we filled x (order filtered)

        path: file, paath for alied fasta 
        od_ident: str

        return: od_ident with the gaps added
        '''
        
        seq_data = find_human_sequence(path)

        # get alied human sequence
        alied_sequence = seq_data['sequence']

        # raise for different length between OMA and uniprot
        protein_id = seq_data['oma_group_id']
        no_gap_sequence_from_oma = seq_data['remove_gap_sequence']
        if len(od_ident) != len(no_gap_sequence_from_oma):
            raise Exception("{}, seq length is different between OMA and uniprot".format(protein_id))

        # make new od_ident
        for index, element in enumerate(alied_sequence):
            if element == "-":
                od_ident = od_ident[:index] + "-" + od_ident[index:]

        # both sides are same in center
        # nen1: 1---1
        for _ in range(2):
            nen1_check = re.finditer("1\-+1", od_ident)
            for i in nen1_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + filed_od * (end - start) + od_ident[end:])

        # nen0: 0---0
        for _ in range(2):
            nen0_check = re.finditer("0\-+0", od_ident)
            for i in nen0_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + filed_od * (end - start) + od_ident[end:])
                    
        # nenx: x---x
        for _ in range(2):
            nenx_check = re.finditer("x\-+x", od_ident)
            for i in nenx_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + filed_od * (end - start) + od_ident[end:])
                    
        # nenz: z---z
        for _ in range(2):
            nenz_check = re.finditer("z\-+z", od_ident)
            for i in nenz_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + filed_od * (end - start) + od_ident[end:])

        # een: ---(0, 1, x, z)
        een_check = re.search("^\-+[0,1,x,z]", od_ident)
        if een_check:
            filed_od = een_check.group()[-1]
            start = een_check.start()
            end = een_check.end()
            od_ident = filed_od * (end - start) + od_ident[end:]
            
        #nee: (0, 1, x, z)---
        nee_check = re.search("[0,1,x,z]\-+$", od_ident)
        if nee_check:
            filed_od = nee_check.group()[0]
            start = nee_check.start()
            end = nee_check.end()
            od_ident = od_ident[:start] + filed_od * (end - start)

        # complicated condition
        # nen01: (0, 1)---(0, 1)
        for _ in range(2):
            nen01_check = re.finditer("[0,1]\-+[0,1]", od_ident)
            od_place = "0" 
            for i in nen01_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + od_place * (end - start) + od_ident[end:])
                    
        # nenxz: (x, z)---(x, z)
        for _ in range(2):
            nenxz_check = re.finditer("[x,z]\-+[x,z]", od_ident)
            od_place = "z" 
            for i in nenxz_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + od_place * (end - start) + od_ident[end:])
                    
        # nen1z: (1, z)---(z, 1)
        for _ in range(2):
            nen1z_check = re.finditer("1\-+z", od_ident)
            od_place = "z" 
            for i in nen1z_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + od_place * (end - start) + od_ident[end:])
        for _ in range(2):
            nen1z_check = re.finditer("z\-+1", od_ident)
            od_place = "z" 
            for i in nen1z_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + od_place * (end - start) + od_ident[end:])
                
        # nen0x: (0, x)---(x, 0)
        for _ in range(2):
            nen0x_check = re.finditer("0\-+x", od_ident)
            od_place = "x" 
            for i in nen0x_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + od_place * (end - start) + od_ident[end:])
        for _ in range(2):
            nen0x_check = re.finditer("x\-+0", od_ident)
            od_place = "x" 
            for i in nen0x_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + od_place * (end - start) + od_ident[end:])
        
        return od_ident
