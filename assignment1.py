#! /usr/bin/env python2
#Attention: Additionally bedtools for Linux has to be installed on the system (sudo apt-get install bedtools)

import mysql.connector
import pysam
import pybedtools

__author__ = 'Priska Lang'


class Assignment1:
    """
        Provides code for working with the bam file AHG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam
            and the gene USH1C.

            This class parses the two files and calculates the following properties using the mysql.connector,
            pysam and pybedtools module:
            - sam header
            - properly paired reads of gene
            - gene reads with indels
            - total average coverage
            - gene average coverage
            - number of mapped reads
            - gene symbols
            - regions of gene
            - numbers of exons
            ! gene average coverage method doesn't work correct!
            For more information see: http://pysam.readthedocs.io/en/latest/api.html, https://daler.github.io/pybedtools/
            and http://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html.
            This code fetches the gene coordinates and provides methods for getting the listed properties above and
            "print_summary" for printing the results.
            """
    def __init__(self):
        """
            The constructor method creates a pysam and a bedtool object from the bam file, fetches the gene USH1C and
            stores it to a dict. Furthermore it stores a list of the properly paired reads and the reads with indels.
            Additionally it counts how often they occur in the bam file. Also the mapped reads of the bam file are
            counted.
            """
        ## Your gene of interest
        self.gene = "USH1C"
        # Read a file in BAM format:
        infileName = "/home/brisi/PycharmProjects/MedizinischeGenomanalysen/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"
        # Create pysam file
        if infileName.endswith(".bam"):
            self.samfile = pysam.AlignmentFile(infileName, "rb")
        # Create pybedtools file
        self.bedtoolfile = pybedtools.BedTool(infileName)
        # Create myGene and fill it with all occurrences of my gene of interest
        self.myGenes = {"name2":[],"name":[],"chrom":[],"start":[],"end":[],"exon":[]}
        self.fetch_gene_coordinates("hg19", "USH1C.TXT")
        # Select the first gene from all occurrences
        self.myGene = {"name2":[],"name":[],"chrom":[],"start":[],"end":[],"exon":[]}
        for key in sorted(self.myGenes.keys()):
            self.myGene[key] = self.myGenes[key][0]
        self.properlyPairedReads = []       # declaration of list for properly paired reads
        self.cProperlyPairedReads = 0       # counter for counting properly paired reads
        self.cMappedReads = 0               # counter for counting mapped reads
        self.geneReadsWithIndels = []       # declaration list of gene reads with indels
        self.cGeneReadsWithIndels = 0       # counter for counting the gene reads with indels
        # determining relevant data in pysam file
        for read in self.samfile.fetch(self.myGene["chrom"], self.myGene["start"], self.myGene["end"]):
            # determining the porperly paired reads
            if read.is_proper_pair:
                self.properlyPairedReads.append(read)
                self.cProperlyPairedReads += 1
            # counting the mapped reads
            if not read.is_unmapped:
                self.cMappedReads += 1
            # determining the reads with indels
            if not read.cigartuples is None and (read.cigartuples[0][0] == 1 or read.cigartuples[0][0] == 2):   # insertion = 1, deletion = 2
                self.geneReadsWithIndels.append(read)
                self.cGeneReadsWithIndels += 1


    def fetch_gene_coordinates(self, genome_reference, file_name):
        """
            Fetches the gene coordinates of genome_reference by using the mysql.connector and saves it as file_name.
            Furthermore it fills the empty dict, created within the constructor method.

            :param genome_reference: the name of the reference genome
            :param file_name: the name of the output file
            """
        print("Connecting to UCSC to fetch data")

        ## Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password',
                                      db=genome_reference)

        ## Get cursor
        cursor = cnx.cursor()

        ## Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]

        ## Build query
        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields)

        ## Execute query
        cursor.execute(query)

        ## Write my gene of interest to file and create it's self.objects
        with open(file_name, "w") as fh:
            for row in cursor:
                if row[0] == self.gene:
                    fh.write(str(row) + "\n")
                    self.myGenes["name2"].append(row[0])
                    self.myGenes["name"].append(row[1])
                    self.myGenes["chrom"].append(row[2])
                    self.myGenes["start"].append(row[3])
                    self.myGenes["end"].append(row[4])
                    self.myGenes["exon"].append(row[6])

        ## Close cursor & connection
        cursor.close()
        cnx.close()

        print("Done fetching data")

    def get_sam_header(self):
        """
            Prints the header of the bam file.
            """
        print("\nheader:")
        print(self.samfile.header["RG"])
        print(self.samfile.header["CO"])
        print(self.samfile.header["HD"])
        # remove comments for printing the whole header object
        #for key in self.samfile.header:
         #   print key
          #  print self.samfile.header[key]
        #todo: embellish print format of sam header

    def get_properly_paired_reads_of_gene(self):
        """
            Returns the properly pared reads of pysam file self.samfile and prints them if you remove the comment signs
            in lines 151-152: #for line in self.properlyPairedReads: and #print line. Additionally it prints the number
            of properly paired reads.

            :return: a list of the properly paired reads of the gene
            """
        # remove comment signs for output on screen
        #for line in self.properlyPairedReads:
        #    print line
        print("\n{} Properly paired reads found.".format(self.cProperlyPairedReads))
        return(self.properlyPairedReads)

    def get_gene_reads_with_indels(self):
        """
            Returns the gene reads with indels from pysam file self.samfile and prints them if you remove the comment
            signs in lines 165-166: #for read in self.geneReadsWithIndels: and #print read. Additionally it prints the
            number of gene reads with indels.

            :return: a list of the gene reads with indels
            """
        # remove comment signs for output on screen
        #for read in self.geneReadsWithIndels:
        #    print read
        print("\n{} gene reads with indels found.".format(self.cGeneReadsWithIndels))
        return(self.geneReadsWithIndels)

    def calculate_average_coverage(self,coverageType):
        """
            Calculates the total average coverage of the bedtool file self.bedtoolfile using the genome_coverage method
            from pybedtools. Afterwards it determines the total average coverage and the gene average coverage and
            prints them if coverageType is "both". If coverageType is "gene", it just prints the genes coverage and
            returns it. Otherwise it prints the total coverage and returns it.

            :param coverageType: the type of coverage that should be printed/returned. Possible values are "both",
            "total" and "gene"

            :return the coverage value if param coverageType is not set to "both"
            """
        self.coverageValues = self.bedtoolfile.genome_coverage(bg=True)
        cbam = 0            # for counting the lines with coverage values in self.bedtoolfile
        cGene = 0           # for counting the lines with coverage values for my gene
        sumValuesBam = 0    # for summing up the coverage values of the self.bedtoolfile
        sumValuesGene = 0   # for summing up the coverage values for my gene
        for line in self.coverageValues:
            cbam += 1
            sumValuesBam += int(line[3])
            # check if the read is within my gene of interest
            if self.myGene["start"] <= line[1] and self.myGene["end"] >= line[2]:
                sumValuesGene += int(line[3])
                cGene += 1
        if coverageType == "both":
            print("\ngene average coverage:")
            print(self.sumValuesGene / float(self.cGene))
            print("\ntotal average coverage:")
            print(float(sumValuesBam) / float(cbam))
        else:
            if coverageType == "gene":
                print("\ngene average coverage:")
                print(self.sumValuesGene / float(self.cGene))
                return(self.sumValuesGene / float(self.cGene))
            else:
                print("\ntotal average coverage:")
                print(float(sumValuesBam) / float(cbam))
                return(float(sumValuesBam) / float(cbam))

    def get_number_mapped_reads(self):
        """
            Prints the number of mapped reads found in the pysam file self.samfile.
            """
        print("\nnumber of mapped reads: ")
        print(self.cMappedReads)

    def get_gene_symbol(self):
        """
            Prints the gene symbol of the first gene in self.myGenes.
            """
        print("\ngene symbol:")
        print(self.myGenes["name2"][0])

    def get_region_of_gene(self):
        """
            Prints the region of gene self.myGene

            :return: a string of gene name, start and end per row, separated with tab
            """
        print("\nregion of gene " + self.myGene["name2"] + ": " + "start: " + str(self.myGene["start"]) + " end: " + \
              str(self.myGene["end"]))
        return(self.myGene["name2"] + "  " + str(self.myGenes["start"]) + " " + str(self.myGenes["end"]) + "   \n")

    def get_number_of_exons(self):
        """
            Prints the list of number of exons of the gene in self.myGene.
            """
        print("\nnumber of exons:")
        print(self.myGene["exon"])

    def print_summary(self):
        """
            print_summary calls all methods above except the constructor method and fetch_gene_coordinates().

                    :Example:

                    For the given input files and not activated printing in self.get_properly_paired_reads_of_gene(),
                    self.get_gene_reads_with_indels() and self.get_region_of_gene() the following is printed on the
                    console:

                :Example:
                Assignment 1
                Connecting to UCSC to fetch data
                Done fetching data

                header:
                [{'LB': '2845856850', 'CN': 'WUGSC', 'DS': 'SRP001294', 'SM': 'HG00096', 'PI': '206', 'ID': 'SRR062634',    'PL': 'ILLUMINA'}, {'LB': '2845856850', 'CN': 'WUGSC', 'DS': 'SRP001294', 'SM': 'HG00096', 'PI': '206', 'ID': 'SRR062635', 'PL': 'ILLUMINA'}, {'LB': '2845856850', 'CN': 'WUGSC', 'DS': 'SRP001294', 'SM': 'HG00096', 'PI': '206', 'ID': 'SRR062641', 'PL': 'ILLUMINA'}]
                ['$known_indels_file(s) = ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.indels.sites.vcf.gz', '$known_indels_file(s) .= ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.low_coverage_vqsr.20101123.indels.sites.vcf.gz', '$known_sites_file(s) = ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.dbsnp.build135.snps.sites.vcf.gz']
                {'SO': 'coordinate', 'VN': '1.0'}

                2382 Properly paired reads found.

                26 gene reads with indels found.

                number of mapped reads:
                2382

                gene symbol:
                USH1C

                region of gene USH1C: start: 17515441 end: 17565963

                number of exons:
                21
                """
        self.get_sam_header()
        self.get_properly_paired_reads_of_gene()
        self.get_gene_reads_with_indels()
        self.calculate_average_coverage("both")
        self.get_number_mapped_reads()
        self.get_gene_symbol()
        self.get_region_of_gene()
        self.get_number_of_exons()

# cue:
if __name__ == '__main__':
    print("Assignment 1")
    assignment1 = Assignment1()
    assignment1.print_summary()