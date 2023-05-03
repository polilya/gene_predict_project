#!/usr/bin/ruby

def printFeature(chromosome, featureName, featureStart, featureEnd, strand, geneName)
  print chromosome,"\tUCSC\t",featureName,"\t",featureStart,"\t",featureEnd,"\t.\t",strand,"\t.\t","gene_id \"#{geneName}\"; transcript_id \"#{geneName}.a\";\n"
end

if (ARGV.size != 1)
  print "Usage: ucsc2gtf.rb <ucsc file>\n"
  exit(0)
end

File.open(ARGV[0], "r") do |ucscFile|
  ucscFile.each_line do |ucscLine|
    lineArray = ucscLine.split
    offset = 0
    if (lineArray.length < 11)
      offset = -1
    end	

    name = lineArray[1 + offset]
    chromosome = lineArray[2 + offset]
    strand = lineArray[3 + offset]
    cdsStart = lineArray[6 + offset].to_i + 1
    cdsEnd = lineArray[7 + offset].to_i
    exonStarts = lineArray[9 + offset]
    exonEnds = lineArray[10 + offset]
    exonStartsArr = exonStarts.split(",")
    exonEndsArr = exonEnds.split(",")
    for i in 0..exonStartsArr.length-1
      start = exonStartsArr[i].to_i + 1
      stop = exonEndsArr[i].to_i
      if (strand == "+")
        if (start < cdsStart)
          if (stop < cdsStart)
            #entire exon is 5' UTR
            printFeature(chromosome, "5UTR", start, stop, strand, name)
          elsif (stop < cdsEnd)
            #partially coding exon on 5' end
            utrEnd = cdsStart - 1
            printFeature(chromosome, "5UTR", start, utrEnd, strand, name)
            printFeature(chromosome, "start_codon", cdsStart, cdsStart+2, strand, name)
            printFeature(chromosome, "CDS", cdsStart, stop, strand, name)
          elsif (stop == cdsEnd)
            #single coding exon
            utrEnd = cdsStart - 1
            printFeature(chromosome, "5UTR", start, utrEnd, strand, name)
            printFeature(chromosome, "start_codon", cdsStart, cdsStart+2, strand, name)
            printFeature(chromosome, "CDS", cdsStart, cdsEnd-3, strand, name)
            printFeature(chromosome, "stop_codon", cdsEnd-2, cdsEnd, strand, name)
          else # stop > cdsEnd
            #single coding exon
            utrEnd = cdsStart - 1
            printFeature(chromosome, "5UTR", start, utrEnd, strand, name)
            printFeature(chromosome, "start_codon", cdsStart, cdsStart+2, strand, name)
            printFeature(chromosome, "CDS", cdsStart, cdsEnd-3, strand, name)
            utrStart = cdsEnd - 2
            printFeature(chromosome, "stop_codon", cdsEnd-2, cdsEnd, strand, name)
            printFeature(chromosome, "3UTR", cdsEnd-2, stop, strand, name)
          end
        elsif (start == cdsStart)
          if (stop < cdsEnd)
            printFeature(chromosome, "start_codon", cdsStart, cdsStart+2, strand, name)
            printFeature(chromosome, "CDS", cdsStart, stop, strand, name)
          elsif (stop == cdsEnd)
            printFeature(chromosome, "start_codon", cdsStart, cdsStart+2, strand, name)
            printFeature(chromosome, "CDS", cdsStart, cdsEnd-3, strand, name)
            printFeature(chromosome, "stop_codon", cdsEnd-2, cdsEnd, strand, name)
          else # stop > cdsEnd
            printFeature(chromosome, "start_codon", cdsStart, cdsStart+2, strand, name)
            printFeature(chromosome, "CDS", cdsStart, cdsEnd-3, strand, name)
            utrStart = cdsEnd - 2
            printFeature(chromosome, "stop_codon", cdsEnd-2, cdsEnd, strand, name)
            printFeature(chromosome, "3UTR", cdsEnd-2, stop, strand, name)
          end
        elsif (start < cdsEnd) # and start > cdsStart
          if (stop < cdsEnd)
            #entire exon is CDS
            printFeature(chromosome, "CDS", start, stop, strand, name)
          elsif (stop == cdsEnd)
            printFeature(chromosome, "CDS", start, cdsEnd-3, strand, name)
            printFeature(chromosome, "stop_codon", cdsEnd-2, cdsEnd, strand, name)
          else # stop > cdsEnd
            printFeature(chromosome, "CDS", start, cdsEnd-3, strand, name)
            printFeature(chromosome, "stop_codon", cdsEnd-2, cdsEnd, strand, name)
            printFeature(chromosome, "3UTR", cdsEnd-2, stop, strand, name)
          end
        else
          #entire exon is 3' UTR
          printFeature(chromosome, "3UTR", start, stop, strand, name)
        end
      else # minus strand
        if (start < cdsStart)
          if (stop < cdsStart)
            #entire exon is 3' UTR
            printFeature(chromosome, "3UTR", start, stop, strand, name)
          elsif (stop < cdsEnd)
            #partially coding exon on 3' end
            utrEnd = cdsStart + 2
            printFeature(chromosome, "3UTR", start, utrEnd, strand, name)
            printFeature(chromosome, "stop_codon", cdsStart, cdsStart+2, strand, name)
            printFeature(chromosome, "CDS", cdsStart+3, stop, strand, name)
          elsif (stop == cdsEnd)
            #single coding exon
            utrEnd = cdsStart + 2
            printFeature(chromosome, "3UTR", start, utrEnd, strand, name)
            printFeature(chromosome, "stop_codon", cdsStart, cdsStart+2, strand, name)
            printFeature(chromosome, "CDS", cdsStart+3, cdsEnd, strand, name)
            utrStart = cdsEnd + 1
            printFeature(chromosome, "start_codon", cdsEnd-2, cdsEnd, strand, name)
          else # stop > cdsEnd
            #single coding exon
            utrEnd = cdsStart + 2
            printFeature(chromosome, "3UTR", start, utrEnd, strand, name)
            printFeature(chromosome, "stop_codon", cdsStart, cdsStart+2, strand, name)
            printFeature(chromosome, "CDS", cdsStart+3, cdsEnd, strand, name)
            utrStart = cdsEnd + 1
            printFeature(chromosome, "start_codon", cdsEnd-2, cdsEnd, strand, name)
            printFeature(chromosome, "5UTR", cdsEnd+1, stop, strand, name)
          end
        elsif (start == cdsStart)
          if (stop < cdsEnd)
            printFeature(chromosome, "stop_codon", cdsStart, cdsStart+2, strand, name)
            printFeature(chromosome, "CDS", cdsStart+3, stop, strand, name)
          elsif (stop == cdsEnd)
            printFeature(chromosome, "stop_codon", cdsStart, cdsStart+2, strand, name)
            printFeature(chromosome, "CDS", cdsStart+3, cdsEnd, strand, name)
            utrStart = cdsEnd + 1
            printFeature(chromosome, "start_codon", cdsEnd-2, cdsEnd, strand, name)
          else # stop > cdsEnd
            printFeature(chromosome, "stop_codon", cdsStart, cdsStart+2, strand, name)
            printFeature(chromosome, "CDS", cdsStart+3, cdsEnd, strand, name)
            utrStart = cdsEnd + 1
            printFeature(chromosome, "start_codon", cdsEnd-2, cdsEnd, strand, name)
            printFeature(chromosome, "5UTR", cdsEnd+1, stop, strand, name)
          end
        elsif (start < cdsEnd) # and start > cdsStart
          if (stop < cdsEnd)
            #entire exon is CDS
            printFeature(chromosome, "CDS", start, stop, strand, name)
          elsif (stop == cdsEnd)
            printFeature(chromosome, "CDS", start, cdsEnd, strand, name)
            utrStart = cdsEnd + 1
            printFeature(chromosome, "start_codon", cdsEnd-2, cdsEnd, strand, name)
          else # stop > cdsEnd
            printFeature(chromosome, "CDS", start, cdsEnd, strand, name)
            utrStart = cdsEnd + 1
            printFeature(chromosome, "start_codon", cdsEnd-2, cdsEnd, strand, name)
            printFeature(chromosome, "5UTR", cdsEnd+1, stop, strand, name)
          end
        else
          #entire exon is 5' UTR
          printFeature(chromosome, "5UTR", start, stop, strand, name)
        end
      end
    end
    print "\n"
  end
end
