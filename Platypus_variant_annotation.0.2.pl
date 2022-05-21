#!/usr/bin/env perl

###########################################################
#### Annotate variants called by platypus
#### sirisha.srao@gmail.com; Wednesday May 18, 2022
####
#### Only works with variants called by Platypus (VCFv4.0)
###########################################################

use strict;
use warnings;
use Getopt::Std;
use JSON;
use Data::Dumper;
use HTTP::Tiny;
use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($INFO);

my $logger = Log::Log4perl->get_logger();
my $http = HTTP::Tiny->new();

my $server = 'http://grch37.rest.ensembl.org';
my $ext = '/vep/human/region';

#######################################################
### Usage
#######################################################

my $usage = qq/

    usage: annotate.pl -i <input_vcf> -o <output_tsv>

/;

die "$usage\n" if (@ARGV != 4);

our ($opt_i, $opt_o, $opt_h);
getopts('i:o:h') || die "$usage\n";

die "$usage\n" if ($opt_h);

######################################################
### Command-line options
######################################################

my $infile = $opt_i;
my $outfile = $opt_o;

open(INVCF, "<$infile") || die "Cannot open $infile $!\n";
open(OUTCSV, ">>$outfile") || die "Cannot open $outfile $!\n";

### print header
print OUTCSV "CHROM,POS,REF,ALT,Depth_of_coverage,variant_read_support,percent_reads_supporting_variant:percent_reads_supporting_reference\n";

my @var_ids = ();


while(<INVCF>) {

    if ($_ !~ /^#/) {                                                       ### work with lines that contain variant information

            my @varfields = split;
            my %info = split(/[;=]/, $varfields[7]);                        ### extract into a hash from the INFO field
            my @alt_alleles = split(/,/, $varfields[4]);                    ### count the number of alternate alleles from the 5th info field.
            my $perc_var;
            my $perc_ref;
            my $var_id;

            #get_var_annotation($varfields[0], $varfields[1], $varfields[4]);

            if (scalar(@alt_alleles) == 1) {
            
                $perc_var = (int($info{'TR'})/int($info{'TC'})) * 100;
                $perc_ref = int($info{'TC'} - $info{'TR'})/int($info{'TC'}) * 100;
                
                print OUTCSV "$varfields[0],$varfields[1],$varfields[3],$varfields[4],$info{'TC'},$info{'TR'}";
                printf OUTCSV ",%3.2f:%3.2f\n", $perc_var, $perc_ref;
                
                $var_id = "$varfields[0]:$varfields[1]-$varfields[1]/$varfields[4]";
                
                push(@var_ids, $var_id);
            
            }
            else {                                                          ### If more than 1 alternate allele exists, record all TR values.

                my @TRs = split(/,/, $info{'TR'});

                for(my $j = 0; $j <= $#TRs; $j++){

                    $perc_var = int($TRs[$j])/int($info{'TC'}) * 100;
                    $perc_ref = int($info{'TC'} - $TRs[$j])/int($info{'TC'}) * 100;

                    print OUTCSV "$varfields[0],$varfields[1],$varfields[3],$alt_alleles[$j],$info{'TC'},$TRs[$j]";
                    printf OUTCSV ",%3.2f:%3.2f\n", $perc_var, $perc_ref;

                    $var_id = "$varfields[0]:$varfields[1]-$varfields[1]/$alt_alleles[$j]";
                    
                    push(@var_ids, $var_id);

                }
                
            }

    }

}

my $var_ids_json_str = encode_json(\@var_ids);
my $jsonfile = "variants.json";

print "$var_ids_json_str\n";

open(VARJSON, ">$jsonfile") || die "Cannot open $jsonfile file to write to $!\n";
print VARJSON $var_ids_json_str;
close VARJSON;


my $response = $http->request('POST', $server.$ext, {
  headers => { 
  	'Content-type' => 'application/json',
  	'Accept' => 'application/json'
  },
  content => '{ "variants" : $var_ids_json_str }'  
  #content => '{ "variants" : ["1:1158631-1158631/G","1:1246004-1246004/G","1:1249187-1249187/A","1:1261824-1261824/C","1:1387667-1387667/G","1:1585597-1585597/G","1:1585642-1585642/T","1:1586752-1586752/C","1:1647686-1647686/C","1:1647722-1647722/TCTAGGATG","1:1647745-1647745/AGCCCTTTT","1:1647778-1647778/G","1:1647893-1647893/CTTTCTT","1:1647968-1647968/CAT","1:1647971-1647971/A","1:1647983-1647983/AGGCTTAT","1:1650787-1650787/CGATGCCTACGTTTC","1:1650807-1650807/C","1:1650832-1650832/G","1:1650845-1650845/A","1:1885055-1885055/T","1:1895177-1895177/T","1:1900106-1900106/TCTC","1:1900232-1900232/C","1:1909843-1909843/T","1:1909854-1909854/A","1:1909868-1909868/A","1:1910001-1910001/G","1:1910012-1910012/C","1:2328714-2328714/G","1:3638674-3638674/T","1:3679775-3679775/T","1:3743109-3743109/G","1:3743132-3743132/C","1:3743391-3743391/T","1:3745787-3745787/C","1:3755675-3755675/C","1:3764216-3764216/T","1:6158562-6158562/G","1:6219287-6219287/TCACA","1:6219287-6219287/T","1:6579521-6579521/C","1:6579607-6579607/T","1:6586070-6586070/T","1:7527996-7527996/C","1:7838113-7838113/C","1:7887493-7887493/C","1:7913445-7913445/T","1:8928125-8928125/T","1:10529194-10529194/C","1:10555257-10555257/T","1:10689924-10689924/T","1:11087524-11087524/A","1:11090916-11090916/A","1:11128001-11128001/G","1:11128654-11128654/A","1:11129848-11129848/G","1:11132217-11132217/A","1:11205058-11205058/T","1:12175729-12175729/T","1:12776218-12776218/C","1:12776344-12776344/T","1:12779560-12779560/C","1:12779618-12779618/C","1:12837634-12837634/A","1:12854530-12854530/G","1:12921332-12921332/CA","1:13173016-13173016/G","1:13940864-13940864/G","1:15447498-15447498/G","1:15687059-15687059/G","1:15812432-15812432/G","1:16259813-16259813/G","1:16373124-16373124/G","1:16380196-16380196/C","1:16382877-16382877/GCCA","1:16383267-16383267/G","1:17012241-17012241/T","1:18554345-18554345/C","1:18554354-18554354/T","1:18554356-18554356/T","1:18702988-18702988/C","1:19411129-19411129/G","1:19413261-19413261/A","1:19415304-19415304/T","1:19474975-19474975/G","1:19519887-19519887/C","1:19565344-19565344/G","1:19595892-19595892/A","1:19596124-19596124/T","1:19596156-19596156/T","1:19611241-19611241/C","1:19948507-19948507/G","1:19950062-19950062/C","1:20216860-20216860/A","1:20233051-20233051/G","1:21167404-21167404/G","1:21754283-21754283/C","1:21795388-21795388/G","1:21797996-21797996/A","1:21799489-21799489/C","1:21801346-21801346/A","1:21805730-21805730/G","1:21805978-21805978/A","1:21935433-21935433/T","1:21940555-21940555/T","1:22013771-22013771/A","1:22183739-22183739/G","1:22216574-22216574/A","1:22216604-22216604/G","1:22217108-22217108/A","1:22329414-22329414/A","1:22408106-22408106/A","1:23399932-23399932/T","1:23418153-23418153/T","1:23419261-23419261/A","1:23419374-23419374/TGGAGGAGCT","1:23419855-23419855/C","1:24394424-24394424/A","1:24394811-24394811/G","1:24400633-24400633/G","1:24485518-24485518/C","1:24658063-24658063/G","1:24927561-24927561/C","1:25362501-25362501/A","1:25362745-25362745/C","1:25362790-25362790/C","1:25570081-25570081/C","1:25678238-25678238/A","1:25780668-25780668/T","1:25780893-25780893/G","1:26084101-26084101/A","1:26152758-26152758/G","1:26162313-26162313/G","1:26310394-26310394/G","1:26357656-26357656/A","1:26383645-26383645/C","1:26524706-26524706/G","1:26609307-26609307/T","1:26646726-26646726/GCATG","1:26655208-26655208/T","1:26664968-26664968/T","1:26771574-26771574/CCATT","1:26879920-26879920/C","1:27427041-27427041/C","1:27682126-27682126/C","1:27995557-27995557/CCCTTACCT","1:27995589-27995589/TGCAC","1:27995605-27995605/T","1:27995873-27995873/A","1:28008847-28008847/G","1:28422913-28422913/T","1:29475341-29475341/T","1:29475394-29475394/G","1:29475648-29475648/G","1:31740706-31740706/G","1:33065947-33065947/C","1:33133968-33133968/C","1:33138516-33138516/C","1:33161212-33161212/C","1:33276424-33276424/A","1:33330432-33330432/G","1:33957152-33957152/G","1:34038214-34038214/C","1:34071525-34071525/T","1:34330067-34330067/C","1:35259961-35259961/G","1:35561510-35561510/A","1:35562965-35562965/A","1:35925861-35925861/A","1:36205166-36205166/G","1:36316571-36316571/C","1:38006359-38006359/G","1:38049552-38049552/G","1:38289383-38289383/C","1:38327982-38327982/T","1:38338795-38338795/G","1:38423853-38423853/C","1:38425724-38425724/C","1:38425895-38425895/C","1:38442547-38442547/A","1:38450367-38450367/G","1:39340862-39340862/G","1:39392668-39392668/T","1:39766041-39766041/G","1:39927507-39927507/A","1:39929240-39929240/C","1:40098328-40098328/C","1:40207169-40207169/G","1:40363054-40363054/C","1:40533266-40533266/G","1:40533287-40533287/G","1:40533347-40533347/G","1:40536074-40536074/CT","1:40735817-40735817/G","1:40879808-40879808/T","1:40881041-40881041/G","1:40998811-40998811/G","1:42049140-42049140/G","1:42898843-42898843/G"] }' 
});

if(!$response->{success}) {
    
    print "Error: $response->{status}\n%$response->{reason}\n";
    while (my ($k, $v) = each %{$response->{headers}}) {    
        for (ref $v eq 'ARRAY' ? @$v : $v) {
            print "$k: $_\n";
        }
    }
    print $response->{content} if length $response->{content};
    
}

if(length $response->{content}) {
  #my $hash = decode_json($response->{content});
  #local $Data::Dumper::Terse = 1;
  #local $Data::Dumper::Indent = 1;
  #print Dumper $hash;
  #print "\n";
}

my $status = $response->{status};


######################################################################
#### subroutines
######################################################################

sub get_var_annotation {

    my ($chr, $coord, $alleles) = @_;

    my @alt_alleles = split(/,/, $alleles);

    foreach my $alt_allele (@alt_alleles) {

        my $ext = "/vep/human/region/$chr:$coord-$coord/$alt_allele?";

        my $response = $http->get($server.$ext, {
            headers => { 'Content-type' => 'application/json' }
        });

        print "Error: ", $response->{status}, "\n" unless $response->{success};

        ### write rate-limiting functionality here.

        print "$response\n";

        if($response->{status} == 429 && exists $response->{headers}->{'retry-after'}) {
             my $retry = $response->{headers}->{'retry-after'};
             sleep($retry);
            
            $response = $http->get($server.$ext, {
            headers => { 'Content-type' => 'application/json' }
            });

        }

        if(length $response->{content}) {
            my $hash = decode_json($response->{content});
            local $Data::Dumper::Terse = 1;
            local $Data::Dumper::Indent = 1;
            print Dumper $hash;
            print "\n\n";
        }
        
    }


}