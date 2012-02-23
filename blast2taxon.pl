#!/usr/bin/perl

#blast2taxon.pl


###############################################################################
#                                                                             #
# Copyright 2012 Anna Friedlander and Monica Gruber.                          #
# (anna.fr@gmail.com and monica.gruber@vuw.ac.nz)                             #
#                                                                             #     
# This program is free software; you can redistribute it and/or modify it     #
# under the same terms as Perl itself (either: the GNU General Public License #
# as published by the Free Software Foundation; or the Artistic License).     #
# See http://dev.perl.org/licenses/ for more information.                     #
#                                                                             #
# Comments and queries are welcome to either or both authors, however we are  #
# unable to provide dedicated technical support.                              #
#                                                                             #
# Please acknowledge the authors if you use the programme in any work         #
# resulting in publication. Please cite it as:                                #
#                                                                             #
#   Friedlander, A.M. and Gruber, M.A.M 2012. blast2taxon.pl v 1.0.           #
#   available from anna.fr@gmail.com or monica.gruber@vuw.ac.nz               #
#                                                                             #
###############################################################################


# blast2taxon.pl is a program to parse BLAST+ output (output format 6 with 
# default fields) and append taxonomic information. 
#
# It contains functions to: 
#    * summarise BLAST+ output by taxon and calculate summary statistics 
#      (number of hits; total, average, and median sequence length; average
#      e-value; average bit-score);
#    * compare BLAST+ output files by number of hits and sequence length at
#      a particular taxon level (or taxon level and GI)
#    * create FASTA files and split data based on whether an input sequence
#      hit a particular taxon.
#
# More detailed information can be found in the help menu by using the command:
#
#      % perl blast2taxon.pl -help

use strict;
use warnings; #comment out to suppress warnings
use Bio::DB::Taxonomy;
use Bio::LITE::Taxonomy::NCBI::Gi2taxid;
use Bio::LITE::Taxonomy::NCBI::Gi2taxid qw/new_dict/;
use Getopt::Long;
use Sort::Naturally;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

#an ordering for taxonomic classification
my @order = qw(species genus subfamily family order class phylum kingdom superkingdom);
my %order_map = map { $order[$_] => $_ } 0 .. $#order;

#command line arguments
my $taxon_info    = '';
my $taxon_summary = '';
my $taxon_overlap = '';
my $to_out        = '';
my $taxon_seq_ol  = '';
my $ts_out        = '';
my $nonunique;
my $level         = '';
my $input_files   = '';
my $sequence_extract;
my $taxon_files   = '';
my $blast_input   = '';
my $name          = '';
my $subseq;
my $help;

my $db;
my $dict;
my $die;

my $nodesfile = '';
my $namesfile = '';
my $dictfile = '';

#parse command line arguments
GetOptions ('taxon_info=s'             => \$taxon_info,
            'taxon_summary=s'          => \$taxon_summary,
            'taxon_overlap=s'          => \$taxon_overlap,
            'to_out=s'                 => \$to_out,
            'taxon_sequence_overlap=s' => \$taxon_seq_ol,
            'ts_out=s'                 => \$ts_out,
            'nonunique'                => \$nonunique,
            'level=s'                  => \$level,
            'nodesfile=s'              => \$nodesfile,
            'namesfile=s'              => \$namesfile,
            'dictfile=s'               => \$dictfile,
            'sequence_extract'         => \$sequence_extract,
            'input_files=s'            => \$input_files,
            'taxon_files=s'            => \$taxon_files,
            'blast_input=s'            => \$blast_input,
            'name=s'                   => \$name,
            'subseq'                   => \$subseq,
            'help|h'                   => \$help);

my @blast_files = split(/\s/,$taxon_info);
my @taxsum      = split(/\s/,$taxon_summary);
my @taxol       = split(/\s/,$taxon_overlap);
my @taxseq      = split(/\s/,$taxon_seq_ol);
my @taxon_files = split(/\s/,$taxon_files);
my @blast_input = split(/\s/,$blast_input);
my @input_files = split(/\s/,$input_files);


#parse command-line arguments and choose subroutine
&choose_sub();


#taxon_info
#fetches taxid and taxonomic hierarchy given GI
#appends to BLAST+ output and prints to file
#(output used in the other four main functions)
sub taxon_info{
    my ($lines, $FNAME, $db, $dict) = @_;

    my @taxidnotfound;
    my @nodenotfound;

    my $OUT = "$FNAME" . "_taxon";
    open (OUT, ">$OUT");

    #print header info
    print OUT "#qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore";
    print OUT "\tgi\ttaxid\tidn\tspecies\tgenus\tsubfam\tfamily\torder\tclass\tphylum\tkingdom\tsuperkingdom\n"; 

    foreach my $line (nsort keys %$lines){
        #get taxid
        my $gi = $lines->{$line}{gi};
        my $tid = $dict->get_taxid($gi);
        if (!$tid){
            push @taxidnotfound, $line;
            delete $lines->{$line};
        }
        else{
            #get node
            my $node = $db->get_taxon(-taxonid => $tid);
            $lines->{$line}{taxid} = $tid;
            #get taxonomic information 
            if($node){
                $lines->{$line}{identifier} = $node->node_name;

                #get lineage information
                while (defined($node)){
                my $rank = $node->rank;
                my $name = $node->node_name;
                $lines->{$line}{"$rank"}= $name;

                #break loop when superkingdom obtained
                if($rank eq "superkingdom"){last;}

                $node = $node -> ancestor;
                }
            }
            else{
               push @nodenotfound, $tid;
            }

        my $idn     = $lines->{$line}{identifier};
        my $sp      = $lines->{$line}{species};
        my $gen     = $lines->{$line}{genus};
        my $sub     = $lines->{$line}{subfamily};
        my $fam     = $lines->{$line}{family};
        my $ord     = $lines->{$line}{order};
        my $cla     = $lines->{$line}{class};
        my $phy     = $lines->{$line}{phylum};
        my $kin     = $lines->{$line}{kingdom};
        my $sup     = $lines->{$line}{superkingdom};

        #print info to file
        printf OUT ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                    ($line,$gi,$tid,$idn,$sp,$gen,$sub,$fam,$ord,$cla,$phy,$kin,$sup));
        }
    }

    close (OUT); 

    #print results for which taxid not found to file
    if(@taxidnotfound){
        $OUT = "$FNAME" . "_taxidnotfound";
        open (OUT, ">$OUT");
        print OUT "# GIs in the NCBI databases may be superceeded with time. To prevent this, ";
        print OUT "ensure you run your BLAST search and this programme against databases ";
        print OUT "from the same date (or as close together in time as is practicable).\n";
        foreach (@taxidnotfound) {print OUT "$_ \n";}
        close (OUT);
    }

    #print results for which taxonomy nodes not found to file
    if(@nodenotfound){
        $OUT = "$FNAME" . "_nodenotfound";
        open (OUT, ">$OUT");
        print OUT "#taxids for which a taxonomy node was not found\n";
        foreach (@nodenotfound) {print OUT "$_ \n";}
        close (OUT);
    }

    return $lines;
}


#taxon_summary
#sorts taxon_info output according to taxon specified by the user
#calculates summary stats (number of hits; total, average, median sequence
#length; average e-value; average bit-score) and prints to file
sub taxon_summary{
    my ($lines,$lvl,$FNAME) = @_;

    #print out headers (column names)
    my $OUT = "$FNAME" . "_summary" . "_$lvl";
    open (OUT, ">$OUT");

    print OUT "#$lvl\thits\ttotal_length\tave_length\tmedian_length\tave_bitscore\tave_evalue";
    foreach my $val (@order){
        if($order_map{$val} > $order_map{$lvl}){
            print OUT "\t$val";
        }
    }
    print OUT "\n";

    #get summary statistics by taxon level specified by user
    my %summary = ();
    foreach my $line (keys %$lines){
        my $taxon = $lines->{$line}{$lvl};
        #deal with special case: taxon is [NONE]
        #(create unique descriptor using upstream information)
        if($taxon eq '[NONE]'){
            my $index = first { $order [$_] eq $lvl } 0 .. $#order;
            my $next = $lvl;
            for(;;){
                last if ($index == @order-1);
                my $next = $lines->{$line}{$order[$index+=1]};
                $taxon = $taxon."-".$next;
                last if (!($next eq '[NONE]'));
            }
        }
        #get upstream tree information (just get once)
        if(!$summary{$taxon}){
            my $str = "";
            foreach my $val (@order){
                if($order_map{$val} > $order_map{$lvl}){
                    $str = $str . "\t" . $lines->{$line}{$val};
                }
            }
        $summary{$taxon}->{tree} = $str;
        @{$summary{$taxon}->{lens}} = ();
        }
        #increment hits and total length, evalue and bitscore info
        $summary{$taxon}->{hits}     += 1;
        $summary{$taxon}->{len}      += $lines->{$line}{len};
        $summary{$taxon}->{evalue}   += $lines->{$line}{evalue};
        $summary{$taxon}->{bitscore} += $lines->{$line}{bitscore};
        push @{$summary{$taxon}->{lens}}, $lines->{$line}{len};
    }

    #print taxon summary data to file
    #NOTE: Length data calculated using ALIGNMENT LENGTH - 
    #ie number of nucleotides that match (not abs(send-sstart))
    for my $key (nsort keys %summary){

        my $value  = $summary{$key};
        my $hits   = $value->{hits};
        my $LEN    = $value->{len};
        my $AVELEN = $LEN/$hits;
        my $MEDLEN = &median(@{$value->{lens}});
        my $AVEBIT = $value->{bitscore}/$hits;
        my $AVEE   = $value->{evalue}/$hits;
        my $TREE   = $value->{tree};

        printf OUT ("%s\t%s\t%s\t%.0f\t%.1f\t%.2f\t%.2e%s\n",
                    ($key,$hits,$LEN,$AVELEN,$MEDLEN,$AVEBIT,$AVEE,$TREE));
    }
    close (OUT); 
}


#taxon_overlap
#sorts taxon_info output according to taxon specified by user
#compares n files by number of hits and total length of hits for
#each taxon
sub taxon_overlap{

    my ($files,$lvl,$n,$OUT,$fnames) = @_;
    my %overlap = ();

    #get data
    my $i = 0;
    foreach my $file (@$files){
        foreach my $k (keys %$file){
            my $taxon = $file->{$k}{$lvl};
            #deal with special case: taxon is [NONE]
            #(create unique descriptor using upstream information)
            if($taxon eq '[NONE]'){
                my $index = first { $order [$_] eq $lvl } 0 .. $#order;
                my $next = $lvl;
                for(;;){
                    last if ($index == @order-1);
                    my $next = $file->{$k}{$order[$index+=1]};
                    $taxon .= "-$next";
                    last if (!($next eq '[NONE]'));
                }
            }
            #get upstream tree information (just get once)
            if(!$overlap{$taxon}){
                my $str = "";
                foreach my $val (@order){
                    if($order_map{$val} > $order_map{$lvl}){
                        $str = $str . "\t" . $file->{$k}{$val};
                    }
                }
            $overlap{$taxon}[$n*2] = $str;
            }
            #add length to total length and increment hits
            $overlap{$taxon}[$i*2]   += $file->{$k}{len};
            $overlap{$taxon}[$i*2+1] += 1;

        }
        $i+=1; 
    }

    #print header information
    open (OUT, ">$OUT");
    open (OUT, ">$OUT");
    print OUT "#$lvl\t";
    for($i=0;$i<$n;$i++){
        printf OUT ("%s length\t%s hits\t", (@$fnames[$i],@$fnames[$i]));
    }
    foreach my $name (@order){
        if($order_map{$name}>$order_map{$lvl}){print OUT "$name\t";}
    }
    print OUT "\n";

    #print data
    foreach my $id (nsort keys %overlap){
        print OUT "$id";
        for($i=0;$i<($n*2);$i++){ 
            my $hits = $overlap{$id}[$i]; if (!$hits){$hits=0;}
            print OUT "\t$hits";
        }   
        printf OUT "%s\n", ($overlap{$id}[$n*2]); 
    }
    close(OUT);

}


#taxon_sequence_overlap
#sorts taxon_info output according to taxon specified by user
#compares n files by number of hits and total length of hits for
#each taxon+GI
#using -nonunique prints only records found in >1 file
sub taxon_sequence_overlap{

    my ($files,$lvl,$n,$OUT,$fnames) = @_;
    my %overlap = ();

    #get data
    my $i=0;
    foreach my $file (@$files){
        foreach my $k (keys %$file){
            my $gi    = $file->{$k}{gi};
            my $taxon = $file->{$k}{$lvl};
            my $id    = $taxon."_".$gi;

            #get classification tree
            if(!$overlap{$id}){
                my $str = "";
                foreach my $val (@order){
                        $str = $str . "\t" . $file->{$k}{$val};
                }
            $overlap{$id}[$n*2] = $str;
            }
            #add length to total length and increment hits
            $overlap{$id}[$i*2]   += $file->{$k}{len};
            $overlap{$id}[$i*2+1] += 1;
        }
        $i+=1; 
    }

    #print header information
    open (OUT, ">$OUT");
    print OUT "#$lvl\_gi\t";
    for($i=0;$i<$n;$i++){
        printf OUT ("%s length\t%s hits\t", (@$fnames[$i],@$fnames[$i]));
    }
    foreach my $name (@order){print OUT "$name\t";}
    print OUT "\n";

    #print data (only print nonunique records, ie found in >1 files)
    if($nonunique){
        foreach my $id (nsort keys %overlap){
            my $filehits = 0;
            for($i=1;$i<($n*2);$i+=2){ 
                if ($overlap{$id}[$i]){$filehits += 1};
            }
            if($filehits > 1){
            print OUT "$id";
                for($i=0;$i<($n*2);$i++){ 
                    my $hits = $overlap{$id}[$i]; if (!$hits){$hits=0;}
                    print OUT "\t$hits";
                }   
                printf OUT "%s\n", ($overlap{$id}[$n*2]); 
            }
        }
    }

    #print data (all data)
    else{
        foreach my $id (nsort keys %overlap){
            print OUT "$id";
            for($i=0;$i<($n*2);$i++){ 
                my $hits = $overlap{$id}[$i]; if (!$hits){$hits=0;}
                print OUT "\t$hits";
            }   
            printf OUT "%s\n", ($overlap{$id}[$n*2]);  
        }
    }
    close(OUT);
}


#sequence_extract
#given taxon_info output and the original input files to BLAST (which
#were used to create the BLAST output used as input to taxon_info),
#and a taxon level and taxon name (e.g. -level superkingdom -name bacteria)
#this function creates two fasta files. One file contains all sequences
#that hit the taxon name specified, the other contains all sequences that 
#hit other taxa. Note that there will likely be overlap in the output files
#(eg an input sequence may hit both bacterial sequences and nonbacterial 
#sequences).
#Using the subseq option prints out only the sequence that was matched
#(using qstart and qend)
sub sequence_extract{

    my ($taxonfiles,$fastafiles,$lvl,$name,$n,$OUTS) = @_;
    print "output file names: \n";

    #separate query sequences based on whether they hit the taxon specified or not
    foreach my $taxon (@$taxonfiles){
        my $fasta = shift @$fastafiles;
        my $OUT   = shift @$OUTS;
        print "$OUT\n$OUT\_nohits\n";
        open (OUT, ">$OUT");
        open (OUT2, ">$OUT\_nohits");
        my %output = ();
        my %nonhits = ();
        #print fasta header and sequence for taxon level and name specified
        foreach my $key (nsort keys %$taxon){
            my $header = $taxon->{$key}{query};
            my $seq    = $fasta->{"$header\n"};
            #print hits
            if(lc($taxon->{$key}{$lvl}) eq lc($name)){
                #if subseq specified, get only sequence that was matched in query
                #and result sequences (using qstart and qend)
                if($subseq){
                    my $qstart = $taxon->{$key}{qstart};
                    my $qend   = $taxon->{$key}{qend};
                    $header = $header."_".$qstart."-".$qend;
                    #strip newlines and get subsequence
                    $seq =~ s/\R*//g;
                    $seq = substr($seq, min($qstart,$qend)-1, abs($qstart-$qend));
                    #add newlines every 50 characters for nice FASTA output
                    $seq =~ s/(.{1,50})/$1\n/gs;
                    }
                #print hit and record in hash (so each hit printed only once)
                if(!$output{$header}){
                    print OUT ">".$name."_".$header."\n".$seq;
                    $output{$header} = $seq;
                }
            }
            #print nonhit and record hash (so each nonhit printed only once)
            else{
                if(!$nonhits{$header}){
                    print OUT2 ">".$name."-nonhit_".$header."\n".$seq;
                    $nonhits{$header} = $seq;
                }
            }
        }
        close(OUT);
        close(OUT2);
    }

}


#help info
if($help){
    print <<HELP;

    USAGE
         perl blast2taxon.pl [-taxon_info filename(s)]
                             [-nodesfile filename]
                             [-namesfile filename]
                             [-dictfile filename]
                             [-taxon_summary use_ti|filename(s)|use_in]
                             [-taxon_overlap use_ti|filename(s)|use_in]
                             [-to_out filename] 
                             [-taxon_sequence_overlap use_ti|filename(s)|use_in] 
                             [-ts_out filename] 
                             [-nonunique]
                             [-level species|genus|subfamily|family|order|class|phylum|kingdom|superkingdom]
                             [-sequence_extract]
                             [-taxon_files filename(s)]
                             [-blast_input filename(s)]
                             [-name taxon_name]
                             [-subseq]
                             [-help]


         -taxon_info appends taxonomic information (species, genus, subfamily, family,
                        order, class) to blast output

                        takes a filename or a space delimited list of filenames in double
                        quotes: "f1 f2 ... fn"

                        files must be in BLAST+ output format 6 (tabular, no 
                        comment lines) with default fields (qseqid sseqid pident length mismatch
                        gapopen qstart qend sstart send evalue bitscore)

                        MUST be used in conjunction with -nodesfile -namesfile and -dictfile

                        OUTPUT: a file <filename_taxon> for each input file which consists of
                        the records in the input file with taxonomic information appended

                        if there are GIs for which taxids are not found, a 
                        <filename_taxidnotfound> containing these GIs (and associated information)
                        will be output

                        if there are taxids for which taxonomy-nodes are not found, a
                        <filename_nodenotfound> containing these taxids will be output


         -nodesfile to specify the NCBI taxonomy nodesfile 

                        takes a filename

                        EG: -nodesfile nodes.dmp

                        get this file from the appropriate taxdump download from 
                        ftp://ftp.ncbi.nih.gov/pub/taxonomy
 
                        use in conjunction with -taxon_info


         -namesfile to specify the NCBI taxonomy namesfile

                        takes a filename

                        EG: -namesfile names.dmp

                        get this file from the appropriate taxdump download from 
                        ftp://ftp.ncbi.nih.gov/pub/taxonomy
 
                        use in conjunction with -taxon_info


         -dictfile to specify the NCBI taxonomy gi->taxid file

                        takes a filename

                        EG: -dictfile gi_taxid_nucl.dmp

                        OR: -dictfile gi_taxid_nucl.bin

                        use the .dmp file the first time you run the programme
                        a .bin file will be created, which you can use subsequently

                        get the appropriate file (nucl or prot) from
                        ftp://ftp.ncbi.nih.gov/pub/taxonomy
 
                        use in conjunction with -taxon_info


         -taxon_summary summarises information by taxon level (hits, total alignment length,
                        average & median alignment length, average bit score, average evalue)

                        takes a filename or a space delimited list of filenames in double
                        quotes: "f1 f2 ... fn"

                        must be in file format as output by -taxon_info

                        alternatively, it can take the argument use_ti if used in 
                        conjunction with -taxon_info

                            EG: -taxon_info <filename> -taxon_summary use_ti

                        if -input_files is used to pass input files, -taxon_summary must be
                        given the argument use_in

                            IE: -taxon_summary use_in -input_files "f1 f2 ... fn"

                        MUST be used in conjunction with -level

                        OUTPUT: a file <filename_taxon_summary_level> for each input file with the
                        following fields: taxon, hits, total alignment length, average
                        alignment length, median alignment length, average bit score, average 
                        evalue, followed by the taxonomic hierarchy upstream from the taxon
                        specified


         -taxon_overlap summarises data by taxon level and file

                        takes a filename or a space delimited list of filenames in double
                        quotes: "f1 f2 ... fn"

                        must be in file format as output by -taxon_info

                        alternatively, it can take the argument use_ti if used in 
                        conjunction with -taxon_info

                            EG: -taxon_info <filename> -taxon_overlap use_ti

                        if -input_files is used to pass input files, -taxon_overlap must be
                        given the argument use_in

                            IE: -taxon_overlap use_in -input_files "f1 f2 ... fn"

                        MUST be used in conjunction with -level and -to_out

                        OUTPUT: one file (<filename> as specified by -to_out) that 
                        summarises the taxon data across all input files, has a taxon 
                        field, followed by the number of hits and total alignment 
                        length for each file, then the taxonomic hierarchy upstream 
                        from the taxon specified


         -to_out use with -taxon_overlap, specifies output file for taxon_overlap

                        takes an output file name


         -taxon_sequence_overlap summarises data by GI, taxon level and file

                        takes a filename or a space delimited list of filenames in double
                        quotes: "f1 f2 ... fn"

                        must be in file format as output by -taxon_info

                        alternatively, it can take the argument use_ti if used in 
                        conjunction with -taxon_info

                            EG: -taxon_info <filename> -taxon_sequence_overlap use_ti

                        if -input_files is used to pass input files, -taxon_overlap must be
                        given the argument use_in

                            IE: -taxon_sequence_overlap use_in -input_files "f1 f2 ... fn"

                        MUST be used in conjunction with -level and -ts_out

                        when used with -nonunique only records found in >1 input file are
                        output

                        OUTPUT: one file (<filename> as specified by -ts_out) that 
                        summarises the taxon data across all input files, fields: 
                        taxon_gi; followed by the number of hits and total alignment 
                        length for each file, then the taxonomic hierarchy upstream 
                        from the taxon specified


         -ts_out use with -taxon_sequence_overlap, specifies out file for taxon_sequence_overlap

                        takes an output file name


         -nonunique when used with -taxon_sequence_overlap only records found in >1 input file
          are output


         -level the taxonomic level to summarise data by

                        takes a keyword from the list: [species|genus|subfamily|family|order|class|phylum|kingdom|superkingdom]


         -sequence_extract takes a list of taxon_files and blast_input files with a taxon level and taxon 
          name (e.g. -level superkingdom -name bacteria) and returns two FASTA files - one containing sequences
          that hit the taxon specified; one that hit other taxa

                        MUST be used in conjunction with -level -name -blast_input and -taxon_files

                        OUTPUT: two files in FASTA format. One file contains all sequences that hit
                        the taxon name specified: <file_seq_name>, the other contains all sequences 
                        that hit other taxa: <file_seq_name_nohits>. Note that there will likely be 
                        overlap in the output files (eg an input sequence may hit both bacterial 
                        sequences and nonbacterial sequences).


         -taxon_files to pass the _taxon files to -sequence_extract

                        takes a file name or a space delimited list of filenames "f1 f2 .. fn"


         -blast_input to pass the blast input files to -sequence_extract (ie the FASTA format
          files used to create the BLAST+ output used as input to taxon_info)

                        takes a file name or a space delimited list of filenames "f1 f2 .. fn"
                        must be in FASTA format

         -name : used to specify a taxon name at the -level specified (e.g. -level superkingdom -name 
          bacteria). For use with sequence_extract.

                        takes a taxon name

         -subseq is an optional argument for sequence_extract. Using subseq means that only the 
          subsequence that matched the specified taxa is printed out to the hits file (using
          qstart and qend)


         -help sends this message to stdout


         ---------NOTE: ERROR MESSAGES--------------------------------------------------------------------------

         An error message like the following:

            substr outside of string at /usr/local/share/perl/5.10.1/Bio/LITE/Taxonomy/NCBI/Gi2taxid.pm line 156.
            Use of uninitialized value in unpack at /usr/local/share/perl/5.10.1/Bio/LITE/Taxonomy/NCBI/Gi2taxid.pm line 156.

         indicates that the taxid was not found for a particular gi (for -taxon_info). This does not affect the execution of
         the programme. The result for which the taxid was not found will be printed out to the _taxid_not_found file. It will
         not be printed out to the _taxon file.

         -------------------------------------------------------------------------------------------------------

HELP
}

#-------PARSE COMMAND-LINE AND CHOOSE SUB-ROUTINE-----------------------------------------------------------------------

sub choose_sub {
    #taxon_info specified
    if(@blast_files){
        if (!$nodesfile | !$namesfile | !$dictfile) {
            print "you must provide a nodesfile, namesfile and gi_taxid file\n\n";
            print "EG: -nodesfile nodes.dmp -namesfile names.dmp -dictfile gi_taxid_nucl.dmp\n\n";
            $die = 1;
        }
        if(!$die){
            #make Bio::DB::Taxonomy database
            print "making taxonomy database ...\n";
            my $db = new Bio::DB::Taxonomy (-source => 'flatfile',
                                            -nodesfile => $nodesfile,
                                            -namesfile => $namesfile);
            print "taxonomy database complete\n";
            #just a little test
            #my $node = $db->get_Taxonomy_Node(-taxonid => 34506);
            #my $node = $db->get_taxon(-taxonid => 34506);
            #if (!$node){print "NODE NOT FOUND\n"; exit;}
    
            my @dictvalues = split(/[.]/,$dictfile);
    
            #make gi_taxid bin file
            if ($dictvalues[1] eq "dmp"){
                print "making gi->taxid bin file ...\n";
                my $dictout = $dictvalues[0] . ".bin";
                new_dict (in =>  $dictfile,
                          out => $dictout);
                print "gi->taxid bin file complete\n";
                $dictfile = $dictout;
            }
    
            #make gi->taxid dictionary
            print "making gi->taxid dictionary ...\n";
            my $dict = Bio::LITE::Taxonomy::NCBI::Gi2taxid->new(dict=>$dictfile);
            print "gi->taxid dictionary complete\n";
        
            #retrieve taxonomic information
            my @files; 
            foreach my $file (@blast_files){
                my $lines = &getGIs($file);
                print "getting taxonomic information for $file...\n";
                $lines = &taxon_info($lines, $file, $db, $dict);
                push @files, $lines;
                printf "got taxonomic information for %s (outfile is %s_taxon)\n",
                       ($file, $file);
    
                #if taxon_summary specified 
                if (@taxsum and $taxsum[0] =~ /^use_ti$/i){
                    my $FNAME = "$file" . "_taxon";
                    &level_check; 
                    print "performing taxon summary for $file...\n"; 
                    &taxon_summary($lines, $level, $FNAME);
                    printf "%s taxon summary for %s complete (outfile is %s_summary_%s)\n",
                           ($level, $file, $FNAME, $level);
                } 
    
            }
            #if taxon_overlap specified 
            if (@taxol and $taxol[0] =~ /^use_ti$/i){
                &level_check; 
                if(!$to_out){ &output_warning('-to_out'); }
                if(!$die){
                    my $n = @files;
                    print "performing taxon overlap...\n";
                    &taxon_overlap(\@files, $level, $n, $to_out, \@blast_files);
                    print "$level taxon overlap complete (outfile is $to_out)\n";
                }
                else{ &die_message; }
            }
            #if taxon_sequence_overlap specified 
            if (@taxseq and $taxseq[0] =~ /^use_ti$/i){
                &level_check; 
                if(!$ts_out){ &output_warning('-ts_out'); }
                if(!$die){
                    my $n = @files;
                    print "performing taxon sequence overlap...\n";
                    &taxon_sequence_overlap(\@files, $level, $n, $ts_out, \@blast_files);
                    print "$level taxon sequence overlap complete (outfile is $ts_out)\n";
                }
                else{ &die_message; }
            }
        }
        else{ &die_message; }
    }
    
    
    #taxon_summary specified
    if (@taxsum and $taxsum[0] !~ /^use_ti$/i){
        #check level and input file(s) are specified
        &level_check; 
        if($taxsum[0] =~ /^use_in$/i and !@input_files){ &input_warning('-taxon_summary'); }
    
        #if all information supplied correctly, proceed; otherwise, die
        if(!$die){
            my @files;
            my $infiles;
            if ($taxsum[0] =~ /^use_in$/i){$infiles = \@input_files;}
            else                          {$infiles = \@taxsum;} 
            foreach my $file (@$infiles){
                my $lines = &getLines($file);
                print "performing taxon summary for $file...\n";
                &taxon_summary($lines, $level, $file);
                printf "%s taxon summary for %s complete (outfile is %s_summary_%s)\n",
                       ($level, $file, $file, $level);
            }
        }
        else{ &die_message; }
    }
    
    #taxon_overlap specified
    if (@taxol and $taxol[0] !~ /^use_ti$/i){
        #check level, input file(s), and output file are specified
        &level_check; 
        if($taxol[0] =~ /^use_in$/i and !@input_files){ &input_warning('-taxon_overlap'); }
        if(!$to_out){ &output_warning('-to_out'); }
    
        #if all information supplied correctly, proceed; otherwise, die
        if(!$die){
            my @files;
            my $infiles;
            if ($taxol[0] =~ /^use_in$/i){$infiles = \@input_files;}
            else                         {$infiles = \@taxol;} 
            foreach my $file (@$infiles){
                my $lines = &getLines($file);
                push @files, $lines;
            }
            my $n = @files;
            print "performing taxon overlap...\n";
            &taxon_overlap(\@files, $level, $n, $to_out, $infiles);
            print "$level taxon overlap complete (outfile is $to_out)\n";
        }
        else{ &die_message; }
    }
    
    #taxon_sequence_overlap specified
    if (@taxseq and $taxseq[0] !~ /^use_ti$/i){
        #check level, input file(s), and output file are specified
        &level_check; 
        if($taxseq[0] =~ /^use_in$/i and !@input_files){ &input_warning('-taxon_sequence_overlap'); }
        if(!$ts_out){ &output_warning('-ts_out'); }
    
        #if all information supplied correctly, proceed; otherwise, die
        if(!$die){
            my @files;
            my $infiles;
            if ($taxseq[0] =~ /^use_in$/i){$infiles = \@input_files;}
            else                          {$infiles = \@taxseq;} 
            foreach my $file (@$infiles){
                my $lines = &getLines($file);
                push @files, $lines;
            }
            my $n = @files;
            print "performing taxon sequence overlap...\n";
            &taxon_sequence_overlap(\@files, $level, $n, $ts_out, $infiles);
            print "$level taxon sequence overlap complete (outfile is $ts_out)\n";
        }
        else{ &die_message; }
    }
    
    #sequence_extract specified
    if ($sequence_extract){
        #check taxon level and taxon name specified; check input files specified
        &level_check; 
        if(!$name){ &specify_warning('a taxon name with -name'); }
        if(!@taxon_files){ &specify_warning('filename_taxon files with -taxon_files');}
        if(!@blast_input){ &specify_warning('blast input files in fasta format with -blast_input');}
        #check taxon_files and blast_input are the same size (ie same number of files)
        if(@taxon_files != @blast_input){
            print "You must specify paired taxon files and blast input files in order\n";
            print "ie: -taxon_files \"file1_taxon .. filen_taxon\" -blast_input \"file1.fa .. filen.fa\"\n\n";
            $die = 1;
        }
        #if all information supplied correctly, proceed; otherwise, die
        if(!$die){
            my @taxonfiles; 
            my @fouts;
            foreach my $file (@taxon_files){
                my $lines = &getLines($file);
                push @taxonfiles, $lines;
                $file =~ s/_taxon/_seq_$name/gi;
                push @fouts, $file;
            }
            my @fastafiles; 
            foreach my $fasta (@blast_input){
                my $lines = &getfasta($fasta);
                push @fastafiles, $lines;
            }
            my $n = @taxonfiles;
            print "performing sequence extract...\n";
            &sequence_extract(\@taxonfiles, \@fastafiles, $level, $name, $n, \@fouts);
            print "sequence extract complete\n";
        }
        else{ &die_message; }
    }

}

#-------END-------------------------------------------------------------------------------------------------


#-------METHODS FOR EXTRACTING DATA FROM FILES--------------------------------------------------------------

#get information from blast ouput file (format 6)
sub getGIs {
    #die if file not found
    open(IN, $_[0]) or die "$_[0] not found. You must specify a blast output file (fmt 6; default fields)\n";
    my %ls = ();
    #get data, put in hash
    while(<IN>){
        chomp;
        my $line = $_;
        my @values = split(/\t/,$_);
        #die if wrong number of fields (will catch some *but not all* incorrect input files)
        if (@values != 12) {die "You must specify a a blast output file (fmt 6; default fields)\n";}
        my $gi = (split(/[|]/,$values[1]))[1];
        my %lineage = ("gi"          , $gi,
                       "query"       , $values[0],
                       "len"         , $values[3],
                       "qstart"      , $values[6],
                       "qend"        , $values[7],
                       "sstart"      , $values[8],
                       "ssend"       , $values[9],
                       "evalue"      , $values[10],
                       "bitscore"    , $values[11],
                       "taxid"       , "[NONE]",
                       "identifier"  , "[NONE]",
                       "species"     , "[NONE]",
                       "genus"       , "[NONE]",
                       "subfamily"   , "[NONE]",
                       "family"      , "[NONE]",
                       "order"       , "[NONE]",
                       "class"       , "[NONE]",
                       "phylum"      , "[NONE]",
                       "kingdom"     , "[NONE]",
                       "superkingdom", "[NONE]");
        $ls{$line} = \%lineage;
    }
    close(IN);
    #return hash ref
    return (\%ls);
}

#get information from fname_taxon file (output by taxon_information)
sub getLines {
    #die if file not found
    open(IN, $_[0]) or die "$_[0] not found. You must specify a file: filename_taxon\n";
    my %ls = ();
    #get data, put in hash
    while(<IN>){
        chomp;
        my $line = $_;
        #ignore comment lines
        if($line !~ /^#/){ 
            my @values = split("\t",$_);
            #die if wrong number of fields (will catch some *but not all* incorrect input files)
            if (@values != 24) {die "You must specify a file: filename_taxon\n";}
            my %lineage = ("query"       , $values[0],
                           "len"         , $values[3],
                           "qstart"      , $values[6],
                           "qend"        , $values[7],
                           "sstart"      , $values[8],
                           "ssend"       , $values[9],
                           "evalue"      , $values[10],
                           "bitscore"    , $values[11],
                           "gi"          , $values[12],
                           "taxid"       , $values[13],
                           "identifier"  , $values[14],
                           "species"     , $values[15],
                           "genus"       , $values[16],
                           "subfamily"   , $values[17],
                           "family"      , $values[18],
                           "order"       , $values[19],
                           "class"       , $values[20],
                           "phylum"      , $values[21],
                           "kingdom"     , $values[22],
                           "superkingdom", $values[23]);
            $ls{$line} = \%lineage;
        }
    }
    close(IN);
    #return hash ref
    return (\%ls);
}

#get information from fasta file
sub getfasta {
    #die if file not found
    open(IN, $_[0]) or die "You must specify a fasta file\n";
    my %ls = ();
    my $header = '';
    my $seq = '';
    #get data, put in hash
    while(<IN>){
        if($_ =~ /^>/){
            if($seq){
                $ls{$header} = $seq;
                $seq = "";
            } 
            $header = substr $_, 1;
        }
        else{
            $seq = $seq . $_;
        }
        $ls{$header} = $seq;
    }
    close(IN);
    #return hash ref
    return (\%ls);
}

#-------END-----------------------------------------------------------------------------------------------------

#-------HELPER METHODS FOR OPTION PARSING-----------------------------------------------------------------------

#if level not specified (correctly), print error message and die
sub level_check{
    if($level !~ /^species$|^genus$|^subfamily$|^family$|^order$|^class$|^phylum$|^kingdom$|^superkingdom$/i){
        print "\nYou must specify a taxon level to sort by with -level\n\n";
        print "USAGE\n  [-level species|genus|subfamily|family|order|class|phylum|kingdom|superkingdom]\n\n";
        $die = 1;
    }
}

#if input file(s) not specified (correctly), print error message and die
sub input_warning{
    printf "You must specify input files (fname_taxon) using %s or -input_files\n\n", ($_[0]);
    printf "USAGE\n  %s \"fname1 .. fnamen\"\n\n", ($_[0]);
    print "   OR\n  -input_files \"fname1 .. fnamen\"\n\n";
    $die = 1;
}

#if output file not specified (correctly), print error message and die
sub output_warning{
    printf "You must specify an output filename with %s\n\n", ($_[0]);
    printf "USAGE\n  [%s filename]\n\n", ($_[0]);
    $die = 1;
}

#if data not specified (correctly), print error message and die
sub specify_warning{
    printf "You must specify %s\n", ($_[0]);
    $die = 1;
}

#print error message and die
sub die_message{
    die "\nCannot proceed due to errors -- programme terminating.\n\n";
}

#-------END HELPER METHODS FOR OPTION PARSING-------------------------------------------------------------------

#helper method to calculate the median value in an array of values
sub median{
    my @a = sort {$a <=> $b} @_;
    return ($a[$#a/2] + $a[@a/2]) / 2;
}