#!/usr/bin/perl

# AmphoraVizu: visualization software for metagenomics analysis tools AMPHORA2 and AmphoraNet

# Copyright (C) 2014 Csaba Kerepesi, Balazs Szalkai, Vince Grolmusz

# AmphoraVizu is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details: <http://www.gnu.org/licenses/>.

# For any other inquiries send an Email to Csaba Kerepesi: kerepesics@gmail.com
# AmphoraVizu webserver: http://pitgroup.org/amphoravizu 

# Options:
# -input: input file (AMPHORA2 output file), default: /dev/stdin
# -min_confidence: minimum confidence, default: 0.1
# -lowest_rank: lowest rank (1=Superkingdom, 7=Species), default: 2
# -min_average: minimum average to display, default: 0.02
# -chart_title: chart title (with optional HTML tags), to override default title
# -chart_type: "column": Column Chart output (default), "pie": Pie chart output with averages
# -chart_width: width of chart area, all CSS values are valid (default: 100%)
# -chart_height: height of chart area (default: 100%)
# -chart_compact: if 1, margins are decreased as much as possible (default: 0)
# -html_inline: if 1, tags like <html>, <head> and <body> are omitted (default: 0)     
# -miniview: 1: miniview on, 0: miniview off (default: 0)      
# -sort: "avg": sorted by averages (default), "abc" sorted by alphabet
# -quantify: "prop": shows the proportions of marker genes (default), "amounts": shows the exact amounts of marker genes

# A typical running command:
# perl AmphoraVizu.pl -input example_output.txt -min_confidence 0.9 -lowest_rank 2 -min_average 0.02 -chart_type column > out.html

my (%avg, %sum, %summarker, %sumphylotype)= ();
my (@phylotypes,@markers);
my $rank="";
my $first=1;
my %args;
my $sum_all_markers=0;

for ($i = 0; $i < @ARGV; $i += 2) {
	$args{$ARGV[$i]} = $ARGV[$i+1];
}
my $infile = ((exists $args{"-input"}) ? $args{"-input"} : "/dev/stdin");
my $min_confidence = ((exists $args{"-min_confidence"}) ? $args{"-min_confidence"} : 0.1);
my $lowest_rank = ((exists $args{"-lowest_rank"}) ? $args{"-lowest_rank"} : 2);
my $min_average = ((exists $args{"-min_average"}) ? $args{"-min_average"} : 0.02);
my $chart_title = ((exists $args{"-chart_title"}) ? $args{"-chart_title"} : "Visualization of \$filename (lowest rank: \$lowest_rank, min. confidence: \$min_confidence, min. average: \$min_average)");
my $chart_type = ((exists $args{"-chart_type"}) ? $args{"-chart_type"} : "column");
my $chart_width = ((exists $args{"-chart_width"}) ? $args{"-chart_width"} : "100%");
my $chart_height = ((exists $args{"-chart_height"}) ? $args{"-chart_height"} : "100%");
my $chart_compact = ((exists $args{"-chart_compact"}) ? $args{"-chart_compact"} : 0);
my $html_inline = ((exists $args{"-html_inline"}) ? $args{"-html_inline"} : 0);     
my $miniview = ((exists $args{"-miniview"}) ? $args{"-miniview"} : 0);              
my $sort = ((exists $args{"-sort"}) ? $args{"-sort"} : "avg");
my $quantify = ((exists $args{"-quantify"}) ? $args{"-quantify"} : "prop");


open (IN, $infile) || die "Can't open $infile";
while (<IN>) {
    if ($first) {
        chop;
        if (! $_ =~ /^Query\s+Marker\s+Superkingdom\s+Phylum\s+Class\s+Order\s+Family\s+Genus\s+Species$/ ) {
        die "Wrong inputfile format. Please use an AMPHORA2 result file as input.\n";
        }
        $first=0;
    } else {
        chop;
        my @inputline;
        my @inputline0 = split /\t/;
        my $j=2;
        my $i=1;
        $inputline[0]=$inputline0[0];
        $inputline[1]=$inputline0[1];
        while ( (exists $inputline0[$j]) && ($j<=1+$lowest_rank) ) {
            if ( $inputline0[$j] =~ /^(.+)\((\d\.\d\d)\)$/ ) {
                if ($2 >= $min_confidence){
                    push @inputline,$1;               
                    push @inputline,$2;
                    $i=$j;
                }
            }
            $j++;
	    }
        if ($i>1) {
            if ($i==2) {
                $rank=""
            } elsif ($i==3) {
                $rank=" Phylum"
            } elsif ($i==4) {
                $rank= " Class"
            } elsif ($i==5) {
                $rank=" Order"
            } elsif ($i==6) {
                $rank=" Family"
            } elsif ($i==7) {
                $rank=" Genus"
            } elsif ($i==8) {
                $rank=""
            }
            my $k=2*$i-2;
            $inputline[$k]=$inputline[$k].$rank;        
            if (!exists $summarker{$inputline[1]}) {
                push @markers,$inputline[1];
            }
            $summarker{$inputline[1]}++;
            if (!exists $sumphylotype{$inputline[$k]}) {
                push @phylotypes, $inputline[$k];
            }
            $sumphylotype{$inputline[$k]}++;
            $sum{$inputline[1]}{$inputline[$k]}++;
        }
    }
}
close IN;

foreach my $keym (@markers) {
    $sum_all_markers+=$summarker{$keym};
}
$nummarkers = @markers;
$num_disp_phylotypes=0;
foreach my $keyp (@phylotypes) {
    foreach my $keym (@markers) {
        if (exists $sum{$keym}{$keyp}) {
            if ($quantify eq "prop") {
                $avg{$keyp}+=$sum{$keym}{$keyp}/$summarker{$keym};
            }elsif ($quantify eq "amounts") {
                $avg{$keyp}+=$sum{$keym}{$keyp};
            }
        }
    }
    if ($quantify eq "prop") {
        $avg{$keyp}=$avg{$keyp}/$nummarkers;
    }elsif ($quantify eq "amounts") {
        $avg{$keyp}=$avg{$keyp}/$sum_all_markers;
    }
    if ($avg{$keyp} > $min_average) {
        $num_disp_phylotypes++;
    }
}
my $calc_width=7*$num_disp_phylotypes*$nummarkers;

if (($chart_type eq "column") && ($miniview) && ($calc_width > 1200)) {     
       printf(qq~Too many columns~);        
}else {                                     
    if (!$html_inline) {            
	    printf("<html><head>");     
	    printf("<meta charset='utf-8'/>");      
	    printf("<title>AmphoraVizu</title>");   
	    printf("</head><body>");        
    }                            
    if($chart_type eq "column") {
        &column_chart();
    }elsif ($chart_type eq "pie") {
        &pie_chart();
    }
    printf(q|<div id="chart_%s" style="width: %s; height: %s;"></div>|, $chart_type, $chart_width, $chart_height);
    if (!$html_inline) {            
        printf("</body>\n");
        printf("</html>\n");
    }                                
}

######################################################################################

sub column_chart {    
    printf(qq~
	    <script type="text/javascript" src="http://www.google.com/jsapi"></script>
	    <script type="text/javascript">
	    google.load('visualization', '1', {packages: ['corechart']});
	    </script>
	    <script type="text/javascript">
	    function drawChart() {
	    var data = google.visualization.arrayToDataTable([
    ~);
    printf("['Phylotypes'"); 
    foreach (@markers) {
	printf(",'%s'",$_);
    }
    printf("]\n");
    if ($sort eq "abc") {
        foreach my $keyp (sort @phylotypes) {
            if ($avg{$keyp} > $min_average) {
                printf(",['$keyp'"); 
                foreach my $keym (@markers) {
                    if (exists $sum{$keym}{$keyp}) { 
                        if ($quantify eq "prop") {
                            printf(",%.2f",$sum{$keym}{$keyp}/$summarker{$keym});
                        }elsif ($quantify eq "amounts") {
                            printf(",%d",$sum{$keym}{$keyp});
                        }
		            } else {
		                printf(",0");
                    }	
	            }
		        printf("]\n");
            }
        }
    } elsif ($sort eq "avg") {
        foreach my $keyp (sort { $avg{$b} <=> $avg{$a} } keys %avg) {
            if ($avg{$keyp} > $min_average) {
                printf(",['$keyp'");
                foreach my $keym (@markers) {
                    if (exists $sum{$keym}{$keyp}) {
                        if ($quantify eq "prop") {
                            printf(",%.2f",$sum{$keym}{$keyp}/$summarker{$keym});
                        }elsif ($quantify eq "amounts") {
                            printf(",%d",$sum{$keym}{$keyp});
                        }
                    } else {
                        printf(",0");
                    }
                }
                printf("]\n");
            }
        }
    }
    printf(qq~  ]);  
	    new google.visualization.ColumnChart(document.getElementById('chart_column')).
	    draw(data, {
    ~);
    if ($miniview == 0) {                   
        $chart_width=1200;
        $chart_height=600;
        if ($chart_width < $calc_width) {
            $chart_width=$calc_width;
        }                                   
    }                                       
    &print_chart_options();
    printf(qq~});
        }
	    google.setOnLoadCallback(drawChart);
        </script>
    ~);
}

sub pie_chart {
    printf(qq~
	    <script type="text/javascript" src="http://www.google.com/jsapi"></script>
	    <script type="text/javascript">
	    google.load('visualization', '1', {packages: ['corechart']});
	    </script>
        <script type="text/javascript">
        function drawChart() {
        var data = google.visualization.arrayToDataTable([
          ['Phylotypes', 'Averages'],
    ~);
    my $j=0;
    foreach my $key (sort { $avg{$b} <=> $avg{$a} } keys %avg) {  
        if ($quantify eq "prop") {
            printf("['%s',%f]",$key,$avg{$key});
        }elsif ($quantify eq "amounts") {
            printf("['%s',%f]",$key,$avg{$key}*$sum_all_markers);
        }
        $j++;
        if ($j != scalar(keys %avg) ) {
            printf(",\n");
        }else {
            printf("\n");
        }
    }
    printf(" ]);\n");
    printf("var options = {");
    &print_chart_options();
    printf(qq~, is3D: true
        };

        var chart = new google.visualization.PieChart(document.getElementById('chart_pie'));
        chart.draw(data, options);
        }
        google.setOnLoadCallback(drawChart);
        </script>
    ~);
}

sub print_chart_options {
    if ($lowest_rank==1) {
        $lowest_rank="Superkingdom"
    } elsif ($lowest_rank==2) {
        $lowest_rank="Phylum"
    } elsif ($lowest_rank==3) {
        $lowest_rank= "Class"
    } elsif ($lowest_rank==4) {
        $lowest_rank="Order"
    } elsif ($lowest_rank==5) {
        $lowest_rank="Family"
    } elsif ($lowest_rank==6) {
        $lowest_rank="Genus"
    } elsif ($lowest_rank==7) {
        $lowest_rank="Species"
    }
    my $title=$chart_title;
    $title =~ s/\$chart_type/$chart_type/eg;
    $title =~ s/\$filename/$infile/eg;
    $title =~ s/\$min_average/$min_average/eg;
    $title =~ s/\$min_confidence/$min_confidence/eg;
    $title =~ s/\$lowest_rank/$lowest_rank/eg;
    $title =~ s/\\/\\\\/g;                          
    $title =~ s/"/\\"/g;                            
    printf(q|title:"%s", width:"%s", height:"%s"|, $title, $chart_width, $chart_height);
    if ($chart_type eq "column") {
        printf(qq~, vAxis: {minValue: 0}~);
    }
    if ($chart_compact) {
        printf(", chartArea: {width:'85%%',height:'80%%'}, legend: {position: 'bottom'}");
    }
}
