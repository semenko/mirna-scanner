#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

my $usage = "STDIN >>> filter_rhybrid.pl <cutoff> >>> STDOUT\n";

die $usage unless @ARGV==1;

my ($maxtar_file,$parameter_file)  = @ARGV;


my %parameters_recent =();
my %parameters =();
my %feature_misc;
my %feature;
my %input_temp =();
my %input;
my %param;
my %input_flot;
my %sites;
my %chisq =();
my %dist  =();
my %param_cutoff =();

my $parameters_recent_ref = 0;
my $dist_ref;
my $source;
my $in = 0;
my $pos_range;
my $parameter;
my $max_tarbase=100;
my $pos_range_string = "1_24|2_7|1_8|1_4|5_8|9_16|9_12|13_16|17_24|17_20|21_24";	
my $seg1_utr_35;
my $seg2_utr_35;
my $seg3_mir_53;
my $seg4_mir_53;
my $parameter_mod='';
my $dist_param;
my $chiprob;
my $total_tarbase;
my $total_rhd;
my $chisq;
my $exist;
my $allcount;
my $newpos;
my $min_chisq =0;
my $param_min_chisq;
my $chisqv;
my $parameters_ref;
my $parameters_temp_ref;
my $utr;
my $mir;
my $energy;
my $pos;
my $utr_len;
my $utr_length;
my $transcript_no;
my $entrez_geneid;
my $pos_orig;
my $alignment_length;
my $parameter_filter;
my $filter_flag;
my $p;
$in     =0;
$source ='rhd';
my $param_cutoff_ref;


my @pos;
my @pos_range;
my @segs;
my @new_pos;
my @parameterList_mod                = (); 
my @chisqlist;
my @seg1_utr_35;



while(<STDIN>)
            {	 
                  $in++;			
	              chomp;
		          ($utr,$utr_length,$mir,undef,$energy,$p,$pos,$seg1_utr_35,$seg2_utr_35,$seg3_mir_53,$seg4_mir_53) = split(/:/);
				  ($entrez_geneid,$transcript_no,undef) =split(/_/,$utr);
				  
				  %feature              = ();
		          %feature_misc         = ();
				  if($source eq 'tarbase')
					{
					        $seg1_utr_35          = reverse($seg1_utr_35);
                            $seg2_utr_35          = reverse($seg2_utr_35);
                            $seg3_mir_53          = reverse($seg3_mir_53); 
                            $seg4_mir_53          = reverse($seg4_mir_53); 
					}
		          scan_alignment($seg1_utr_35, $seg2_utr_35, $seg3_mir_53, $seg4_mir_53,\%feature,$in,$source);			    
		          my $feature_misc_ref  = scan_alignment4misc_params($seg1_utr_35,$seg2_utr_35,$seg3_mir_53,$seg4_mir_53,\%feature,$in,$source);
                  my %feature_misc      = %$feature_misc_ref;
		          $parameters_recent_ref= pass_filter(\%feature,\%feature_misc,\%param,$source,$pos_range_string,$in);# returned paparmeter(s);	
				  %parameters_recent    = %$parameters_recent_ref;
				  $utr_len              = get_utrlength($seg1_utr_35, $seg2_utr_35 );
				  $param_cutoff_ref     = get_param_cutoff();
				  %param_cutoff         = %$param_cutoff_ref;
				  
                  $pos_orig             = $utr_length-$pos-$utr_len+2;
				  
				  $filter_flag          = 0;
				  
				  foreach $parameter_filter (sort keys %param_cutoff)
                   {
	                    if($parameters_recent{$source}{$in}{$parameter_filter}==-1)
	                      {
	                          $parameters_recent{$source}{$in}{$parameter_filter} =0;
	                      }
	                    else
	                      {
	                      }
						if(($parameter_filter eq 'matchcount1_8') or ($parameter_filter eq 'matchcount1_12') or ($parameter_filter eq 'matchcount1_8overmatchcount17_24'))
						{
                          if($parameters_recent{$source}{$in}{$parameter_filter}<$param_cutoff{$parameter_filter})
	                       {		
		                     $filter_flag++;
		                   }
						else
						  {						     
						  }
						}
					   else
					   {
					     if($parameters_recent{$source}{$in}{$parameter_filter}>$param_cutoff{$parameter_filter})
	                       {		
		                     $filter_flag++;
		                   }
						 else
						   {
						   }
					  }
                   }
                  if($filter_flag==0)
                   {
                         $seg1_utr_35 = reverse($seg1_utr_35);
	                     $seg2_utr_35 = reverse($seg2_utr_35);
	                     $seg3_mir_53 = reverse($seg3_mir_53);
	                     $seg4_mir_53 = reverse($seg4_mir_53);
						 print "$entrez_geneid\t$transcript_no\t$mir\t$energy\t$pos_orig\t$seg1_utr_35\t$seg2_utr_35\t$seg3_mir_53\t$seg4_mir_53\n";
				   }
                  else
                  {
                  }
				%parameters_recent    = ();
             }


      
sub get_param_cutoff
{
 my %param_cutoff;
 open(CUT,$maxtar_file) || die ("Could not open $maxtar_file");
 while(<CUT>)
 {
    if($_=~m/^#/)
	{
	    #Dont do anything
	}
    else
	{
	   my($param_tar,$maxtar)    = split(/\t/);
	   $param_cutoff{$param_tar} = $maxtar;
	}
 }
 close(CUT);
 return \%param_cutoff;
}	 
	 
	 
  
sub get_utrlength
{
	my ($seg1_utr_35, $seg2_utr_35 ) = @_;
	my $utr_len;
	for (my $i=0; $i<length($seg1_utr_35); $i++)
	{
		if ( substr( $seg1_utr_35, $i, 1 ) ne ' ' || substr( $seg2_utr_35, $i, 1 ) ne ' ' )
		{
			$utr_len++;
		}
	}
	return $utr_len;
}


sub pass_filter
{    
    my ($feature_ref, $feature_misc_ref,$param_ref,$source,$pos_range_string,$in) = @_;
	
	my $params;
	my @parameterModified     = ();
	my @parameterListModified = ();
	my %temphash              = ();
	my %feature_misc          = %$feature_misc_ref; 
	my %parameters_recent            = ();
		
        @pos_range                           = split(/\|/,$pos_range_string);
	    my $gaps1_8overgaps17_24;
	    my $matchcount1_8overmatchcount9_16;
	    my $matchcount1_12overmatchcount13_24;
        my $mismatchbasecount1_12overmismatchbasecount13_24;
        my $matchcount1_24overmismatchbasecount1_24;
        my $matchcount1_12overmismatchbasecount1_12;
        my $matchcount13_24overmismatchbasecount13_24;
        my $wcpairsoverwobble;
        my $terminalbasesoverinnerbases;
        my $matchcount1_8overmatchcount17_24;
        my $matchcount9_16overmatchcount17_24;
        my $mismatchbasecount1_8overmismatchbasecount9_16;
        my $mismatchbasecount1_8overmismatchbasecount17_24;
        my $mismatchbasecount9_16overmismatchbasecount17_24;
        my $gaps1_8overgaps9_16;
        my $gaps9_16overgaps17_24;
        my $loop1_8overloop9_16;
        my $loop1_8overloop17_24;
        my $loop9_16overloop17_24;
        my $mbulge1_8overmbulge9_16;
        my $mbulge1_8overmbulge17_24;
        my $mbulge9_16overmbulge17_24;
        my $ubulge1_8overubulge9_16;
        my $ubulge1_8overubulge17_24 ;
        my $ubulge9_16overubulge17_24;
        my $mismatch1_8overmismatch17_24;
        my $rangecount=0;
    foreach my $pos_range(@pos_range)
     {
        my ($posBegin,$posEnd)            = split(/_/,$pos_range);
       
	    my $rangecount++;      
          
        my @loop =();
	
	   my @mbulge=();
	
	   my @ubulge=();
	
	   my @loop_4mismatch=();
	
	   my @mbulge_4mismatch=();
	
	   my @ubulge_4mismatch=();

	   my @assymetric_loops=();
	   
	   my @looplength_mir=();

	   my @looplength_utr=();
	   
	   my  @looplens=(); #This holds mir and utr lens in loop
	   
	
	   
	   my $gu = 0;
	
	   my $gu_4all = 0;
	   
	   my $term_mismatch;

	   my $term_mismatch_mir35;
	
	   my $oversize_loop = 0;
	
	   my $oversize_mbulge = 0;
	
	   my $oversize_ubulge = 0;

	   my $matchseglen;

	   my $matchcount;

	   my $mismatch; # ubulge or mbulge or loop
	
	   my $all; # ubulge or mbulge or loop or term mismatch (>2) or GU
	
	  
	
	
	foreach my $gu_ref (@{$feature_ref->{gu}})
	{
		if( $gu_ref->{mirpos}+1 >= $posBegin && $gu_ref->{mirpos}+1 <= $posEnd)
		{
			$gu++;
		}
	}
	
	foreach my $gu_ref (@{$feature_ref->{gu}})
	{
		if( $gu_ref->{mirpos}+1 >= $posBegin && $gu_ref->{mirpos}+1 <= $posEnd )
		{
			$gu_4all++;
		}
	}

	########################################################################################################################
	#   correct only if the range starts at 1 otherwise algorithm needs to be modified
	########################################################################################################################

	foreach my $loop_ref (@{$feature_ref->{loop}})
	{
		if( $loop_ref->{miropen}+1 >= $posBegin && $loop_ref->{miropen}+1 < $posEnd )
		 {           
                push @looplength_mir,$loop_ref->{mirclose} - $loop_ref->{miropen} - 1 ;
                       
		    push @looplength_utr,$loop_ref->{utrclose} - $loop_ref->{utropen} - 1;
		    
		
		    if( $loop_ref->{mirclose} - $loop_ref->{miropen} - 1 != $loop_ref->{utrclose} - $loop_ref->{utropen} - 1 )
			{
                push @assymetric_loops,$loop_ref->{mirclose} - $loop_ref->{miropen} - 1;
			}			    
            
			if( $loop_ref->{mirclose} - $loop_ref->{miropen} - 1 > $loop_ref->{utrclose} - $loop_ref->{utropen} - 1 )
			{
				push @loop, $loop_ref->{mirclose} - $loop_ref->{miropen} - 1;
			}
			else
			{
				push @loop, $loop_ref->{utrclose} - $loop_ref->{utropen} - 1;			
			}
			
		}
	}
	foreach my $mbulge_ref (@{$feature_ref->{mbulge}})
	{
		if( $mbulge_ref->{miropen}+1 >= $posBegin && $mbulge_ref->{miropen}+1 < $posEnd )
		{
			push @mbulge, $mbulge_ref->{mirclose} - $mbulge_ref->{miropen} - 1;
	    }
    }
	
	foreach my $ubulge_ref (@{$feature_ref->{ubulge}})
	 {
		if($ubulge_ref->{miropen}+1 >= $posBegin && $ubulge_ref->{miropen}+1 < $posEnd )
		  {
			push @ubulge, $ubulge_ref->{utrclose} - $ubulge_ref->{utropen} - 1;

		  }
       }

	# for calculating mismatch

	foreach my $loop_ref (@{$feature_ref->{loop}})
	  {
		if( $loop_ref->{miropen}+1 >= $posBegin && $loop_ref->{miropen}+1 < $posEnd )
		 {
			if( $loop_ref->{mirclose} - $loop_ref->{miropen} - 1 > $loop_ref->{utrclose} - $loop_ref->{utropen} - 1 )
			 {
				push @loop_4mismatch, $loop_ref->{mirclose} - $loop_ref->{miropen} - 1;
			 }
			else
			 {
				push @loop_4mismatch, $loop_ref->{utrclose} - $loop_ref->{utropen} - 1;
			 }
		 }
	 }
	
	foreach my $mbulge_ref (@{$feature_ref->{mbulge}})
	 {
		if( $mbulge_ref->{miropen}+1 >= $posBegin && $mbulge_ref->{miropen}+1 < $posEnd )
		 {
			push @mbulge_4mismatch, $mbulge_ref->{mirclose} - $mbulge_ref->{miropen} - 1;
		 }
	 }
	
	foreach my $ubulge_ref (@{$feature_ref->{ubulge}})
	 {
		if( $ubulge_ref->{miropen}+1 >= $posBegin && $ubulge_ref->{miropen}+1 < $posEnd )
		{
			push @ubulge_4mismatch, $ubulge_ref->{utrclose} - $ubulge_ref->{utropen} - 1;
		}
	 }

	########################################################################################################################

	$mismatch = scalar(@loop_4mismatch) + scalar(@mbulge_4mismatch) + scalar(@ubulge_4mismatch);
	
	$term_mismatch               = $feature_ref->{term_mismatch};

	$term_mismatch_mir35         = $feature_ref->{term_mismatch_mir35};   
	$matchseglen                 = $feature_ref->{matchseglen};
	
    my $assymetricloopcount      = scalar(@assymetric_loops);
	
	my $term_mismatch_penalty    = $term_mismatch > 1? $term_mismatch-1: 0;
	
	$all                         = $mismatch + $gu_4all + $term_mismatch_penalty; 
		
    $matchcount                  = $feature_ref->{matchcount};   
        
    my $scalar_mbulge            = scalar(@mbulge);
         
    my $scalar_ubulge            = scalar(@ubulge);
        
    my $scalar_loop              = scalar(@loop);

  
    if(scalar(@loop)              ==0) {push @loop,'-1';              }          
    if(scalar(@ubulge)            ==0) {push @ubulge,'-1';            }
    if(scalar(@mbulge)            ==0) {push @mbulge,'-1';            }
    if(scalar(@looplength_mir)    ==0) {push @looplength_mir,'-1';    }
    if(scalar(@looplength_utr)    ==0) {push @looplength_utr,'-1';    }
  

   @loop                      = sort {$b <=> $a} @loop;
   @ubulge                    = sort {$b <=> $a} @ubulge;   	
   @mbulge                    = sort {$b <=> $a} @mbulge;
   @looplength_mir            = sort {$b <=> $a} @looplength_mir;
   @looplength_utr            = sort {$b <=> $a} @looplength_utr;


   my $max_looplength         = shift(@loop);
   my $max_ubulge             = shift(@ubulge);
   my $max_mbulge             = shift(@mbulge);
   my $max_looplength_mir     = shift(@looplength_mir);
   my $max_looplength_utr     = shift(@looplength_utr);


   my @parameters             = ($scalar_loop,$scalar_mbulge,$scalar_ubulge,,$gu,$mismatch,$all,$term_mismatch,$max_looplength,$max_ubulge,
                                 $max_mbulge,$max_looplength_mir,$max_looplength_utr,$term_mismatch_mir35,$matchseglen,
  								 $assymetricloopcount);
   my @parameterList          = ('loop','mbulge','ubulge','gu','mismatch','all','term_mismatch','max_looplength','max_ubulge','max_mbulge',
                                 'max_looplength_mir','max_looplength_utr','term_mismatch_mir35','matchseglen',
								 'assymetricloopcount');
   for(my $i=0;$i<scalar(@parameters);$i++)
    {  
	    if(($rangecount>1) &&($parameterListModified[$i] eq 'term_mismatch' || $parameterListModified[$i] eq 'term_mismatch_mir35' || $parameterListModified[$i] eq 'matchseglen'))
            {
			   next;
	        }
		my $parameterListmodified                           = $parameterList[$i].$posBegin.'_'.$posEnd;
		push @parameterModified,$parameters[$i];
		push @parameterListModified,$parameterListmodified;
    }		
  }
  
  


  foreach  my $param_misc( keys %feature_misc)           
   {
	  
      if($param_misc eq 'mir_length')
	    {
	       next;
	    }
     else
	    {
           push @parameterModified, $feature_misc{$param_misc};
	       push @parameterListModified,$param_misc;
	    }
   }
 
 


for(my $u =0;$u<scalar(@parameterModified);$u++)
 {
    $temphash{$parameterListModified[$u]} = $parameterModified[$u];
 }




  if($temphash{'matchcount13_24'}<=0)
	{
	    $matchcount1_12overmatchcount13_24 = -2;
	}
    else
	{  
	   $matchcount1_12overmatchcount13_24                = $temphash{'matchcount1_12'}/$temphash{'matchcount13_24'};
	}
	
    push @parameterModified,$matchcount1_12overmatchcount13_24;
    push @parameterListModified,'matchcount1_12overmatchcount13_24';

  if($temphash{'mismatchbasecount13_24'}<=0) 
   {
     $mismatchbasecount1_12overmismatchbasecount13_24 =-2;
   }
  else
   {
     $mismatchbasecount1_12overmismatchbasecount13_24      = $temphash{'mismatchbasecount1_12'}/$temphash{'mismatchbasecount13_24'};
   }
     push @parameterModified,$mismatchbasecount1_12overmismatchbasecount13_24;
     push @parameterListModified,'mismatchbasecount1_12overmismatchbasecount13_24;';
	
 
  if($temphash{'mismatchbasecount1_24'}<=0)
   {
     $matchcount1_24overmismatchbasecount1_24 =-2;
   }
 else
   {
     $matchcount1_24overmismatchbasecount1_24         = $temphash{'matchcount1_24'}/$temphash{'mismatchbasecount1_24'};
   }
push @parameterModified,$matchcount1_24overmismatchbasecount1_24;
push @parameterListModified,'matchcount1_24overmismatchbasecount1_24';

if($temphash{'mismatchbasecount1_12'}<=0)
{
   $matchcount1_12overmismatchbasecount1_12=-2;
}
else
{
   $matchcount1_12overmismatchbasecount1_12         = $temphash{'matchcount1_12'}/$temphash{'mismatchbasecount1_12'};
}
push @parameterModified,$matchcount1_12overmismatchbasecount1_12;
push @parameterListModified,'matchcount1_12overmismatchbasecount1_12';
if
($temphash{'mismatchbasecount13_24'}<=0)
{
  $matchcount13_24overmismatchbasecount13_24 =-2;
}
else
{
 $matchcount13_24overmismatchbasecount13_24         = $temphash{'matchcount13_24'}/$temphash{'mismatchbasecount13_24'};
}
push @parameterModified,$matchcount13_24overmismatchbasecount13_24;
push @parameterListModified,'matchcount13_24overmismatchbasecount13_24';

if($temphash{'gu1_24'}<=0)
{
 $wcpairsoverwobble=-2;
}
else
{
 $wcpairsoverwobble = $feature_misc{'mir_length'}/$temphash{'gu1_24'};
}
push @parameterModified,$wcpairsoverwobble;
push @parameterListModified,'wcpairsoverwobble';
                            
if(($feature_misc{'mir_length'}-$feature_ref->{'term_mismatch_mir35'}-$feature_ref->{'term_mismatch'})<=0)
{
  $terminalbasesoverinnerbases =-2;
}
else
{
  $terminalbasesoverinnerbases = ($feature_ref->{'term_mismatch_mir35'}+$feature_ref->{'term_mismatch'})/($feature_misc{'mir_length'}-$feature_ref->{'term_mismatch_mir35'}-$feature_ref->{'term_mismatch'});
}
push @parameterModified,$terminalbasesoverinnerbases;
push @parameterListModified,'terminalbasesoverinnerbases';

if($temphash{'matchcount9_16'}<=0)
{
  $matchcount1_8overmatchcount9_16 = -2;
}
else
{
  $matchcount1_8overmatchcount9_16                 = $temphash{'matchcount1_8'}/$temphash{'matchcount9_16'};
}
push @parameterModified,$matchcount1_8overmatchcount9_16;
push @parameterListModified,'matchcount1_8overmatchcount9_16';

if($temphash{'matchcount17_24'}<=0)
{
  $matchcount1_8overmatchcount17_24 = $temphash{'matchcount1_8'};
}
else
{
  $matchcount1_8overmatchcount17_24                = $temphash{'matchcount1_8'}/$temphash{'matchcount17_24'};
}
push @parameterModified,$matchcount1_8overmatchcount17_24;
push @parameterListModified,'matchcount1_8overmatchcount17_24';

if($temphash{'matchcount17_24'}<=0)
{
  $matchcount9_16overmatchcount17_24=-2;
}
else
{
  $matchcount9_16overmatchcount17_24                = $temphash{'matchcount9_16'}/$temphash{'matchcount17_24'};
}

push @parameterModified,$matchcount9_16overmatchcount17_24;
push @parameterListModified,'matchcount9_16overmatchcount17_24';

if($temphash{'mismatchbasecount9_16'}<=0)
{
  $mismatchbasecount1_8overmismatchbasecount9_16=-2;
}
else
{
  $mismatchbasecount1_8overmismatchbasecount9_16                = $temphash{'mismatchbasecount1_8'}/$temphash{'mismatchbasecount9_16'};
}
push @parameterModified,$mismatchbasecount1_8overmismatchbasecount9_16;
push @parameterListModified,'mismatchbasecount1_8overmismatchbasecount9_16';

if($temphash{'mismatchbasecount17_24'}<=0)
{
 $mismatchbasecount1_8overmismatchbasecount17_24=-2;
}
else
{
 $mismatchbasecount1_8overmismatchbasecount17_24                = $temphash{'mismatchbasecount1_8'}/$temphash{'mismatchbasecount17_24'};
}
push @parameterModified,$mismatchbasecount1_8overmismatchbasecount17_24;
push @parameterListModified,'mismatchbasecount1_8overmismatchbasecount17_24';

if($temphash{'mismatchbasecount17_24'}<=0)
{
 $mismatchbasecount9_16overmismatchbasecount17_24=-2;
}
else
{ 
 $mismatchbasecount9_16overmismatchbasecount17_24                = $temphash{'mismatchbasecount9_16'}/$temphash{'mismatchbasecount17_24'};
}
push @parameterModified,$mismatchbasecount9_16overmismatchbasecount17_24;
push @parameterListModified,'mismatchbasecount9_16overmismatchbasecount17_24';

if($temphash{'gap9_16'}<=0)
{
 $gaps1_8overgaps9_16=-2;
}
else
{ 
 $gaps1_8overgaps9_16                = $temphash{'gap1_8'}/$temphash{'gap9_16'};
}  

push @parameterModified,$gaps1_8overgaps9_16;
push @parameterListModified,'gaps1_8overgaps9_16';
                       
if($temphash{'gap17_24'}<=0)
 {
  $gaps1_8overgaps17_24=-2;
 }
else
{ 
 $gaps1_8overgaps17_24                = $temphash{'gap1_8'}/$temphash{'gap17_24'};
}

push @parameterModified,$gaps1_8overgaps17_24;
push @parameterListModified,'gaps1_8overgaps17_24';

if($temphash{'gap17_24'}<=0)
{
 $gaps9_16overgaps17_24=-2;
}
else
{
 $gaps9_16overgaps17_24                = $temphash{'gap9_16'}/$temphash{'gap17_24'};
}
push @parameterModified,$gaps9_16overgaps17_24;
push @parameterListModified,'gaps9_16overgaps17_24';

if($temphash{'loop9_16'}<=0)
{
 $loop1_8overloop9_16=-2;
}
else
{
 $loop1_8overloop9_16                = $temphash{'loop1_8'}/$temphash{'loop9_16'};
}

push @parameterModified,$loop1_8overloop9_16;
push @parameterListModified,'loop1_8overloop9_16';
if($temphash{'loop17_24'}<=0)
{
 $loop1_8overloop17_24=-2;
}
else
{
 $loop1_8overloop17_24                = $temphash{'loop1_8'}/$temphash{'loop17_24'};
}

push @parameterModified,$loop1_8overloop17_24;
push @parameterListModified,'loop1_8overloop17_24';
	

if($temphash{'loop17_24'}<=0)
{
 $loop9_16overloop17_24=-2;
}
else
{
 $loop9_16overloop17_24                = $temphash{'loop9_16'}/$temphash{'loop17_24'};
}

push @parameterModified,$loop9_16overloop17_24;
push @parameterListModified,'loop9_16overloop17_24';
	
if($temphash{'mbulge9_16'}<=0)
{
 $mbulge1_8overmbulge9_16=-2;
}
else
{
 $mbulge1_8overmbulge9_16                = $temphash{'mbulge1_8'}/$temphash{'mbulge9_16'};
}

push @parameterModified,$mbulge1_8overmbulge9_16;
push @parameterListModified,'mbulge1_8overmbulge9_16';

if($temphash{'mbulge17_24'}<=0)
{
 $mbulge1_8overmbulge17_24=-2;
}
else
{
 $mbulge1_8overmbulge17_24               = $temphash{'mbulge1_8'}/$temphash{'mbulge17_24'};
}

push @parameterModified,$mbulge1_8overmbulge17_24;
push @parameterListModified,'mbulge1_8overmbulge17_24';
	
if($temphash{'mbulge17_24'}<=0)
{
 $mbulge9_16overmbulge17_24=-2;
}
else
{
 $mbulge9_16overmbulge17_24                = $temphash{'mbulge9_16'}/$temphash{'mbulge17_24'};
}

push @parameterModified,$mbulge9_16overmbulge17_24;
push @parameterListModified,'mbulge9_16overmbulge17_24';
	
if($temphash{'ubulge9_16'}<=0)
{
$ubulge1_8overubulge9_16=-2;
}
else
{
$ubulge1_8overubulge9_16                = $temphash{'ubulge1_8'}/$temphash{'ubulge9_16'};
}

push @parameterModified,$ubulge1_8overubulge9_16;
push @parameterListModified,'ubulge1_8overubulge9_16';

if($temphash{'matchcount17_24'}<=0)
{
$ubulge1_8overubulge17_24 =-2;
}
else
{ 
$ubulge1_8overubulge17_24                 = $temphash{'matchcount1_8'}/$temphash{'matchcount17_24'};
}

push @parameterModified,$ubulge1_8overubulge17_24;
push @parameterListModified,'ubulge1_8overubulge17_24';
	

if($temphash{'ubulge17_24'}<=0)
{
$ubulge9_16overubulge17_24=-2;
}
else
{ 
$ubulge9_16overubulge17_24                = $temphash{'ubulge9_16'}/$temphash{'ubulge17_24'};
}

push @parameterModified,$ubulge9_16overubulge17_24;
push @parameterListModified,'ubulge9_16overubulge17_24';

if($temphash{'mismatch17_24'}<=0)
{
$mismatch1_8overmismatch17_24=-2;
}
else
{ 
$mismatch1_8overmismatch17_24                = $temphash{'mismatch1_8'}/$temphash{'mismatch17_24'};
} 

push @parameterModified,$mismatch1_8overmismatch17_24;
push @parameterListModified,'mismatch1_8overmismatch17_24';
%temphash   =();

 for(my $i =0;$i<scalar(@parameterListModified);$i++)
  {
      if($parameterListModified[$i]=~m/(\S+)over(\S+)/)	
	   {
         #Do nothing
	   }
     elsif($parameterListModified[$i]=~m/^term_mismatch/)
           {
            if(($parameterListModified[$i] eq 'term_mismatch1_24') || ($parameterListModified[$i] eq 'term_mismatch_mir351_24'))
                {
                    $parameters_recent{$source}{$in}{$parameterListModified[$i]} = $parameterModified[$i];
                }
	    else
		{
				   next;
		}
           }
    elsif($parameterListModified[$i]=~m/^matchseglen/)
          {
            if($parameterListModified[$i] eq 'matchseglen1_24')
                {
                    $parameters_recent{$source}{$in}{$parameterListModified[$i]} = $parameterModified[$i];
                }
	   else
		{
	            next;
	        }
          }
		
      $parameters_recent{$source}{$in}{$parameterListModified[$i]} = $parameterModified[$i];
  }  
   
  return \%parameters_recent;
} 



sub get_param
{
 	my ($parameter_file, $param_ref) = @_;
 	
 	open (PARAM,$parameter_file) || die ("Could not open $parameter_file!"); 
	
	while(<PARAM>)
	{
		if($_=~m/gupair allowed\t(\d+)/)
		{
			$param_ref->{gu_count} = $1; 
		}
		elsif($_=~m/gupair range\t(\d+)__(\d+)/)
		{
			( $param_ref->{gu_begin}, $param_ref->{gu_end} ) = ($1,$2);
		}
		elsif($_=~m/mir-bulge allowed\t(\d+)/)
		{ 
		$param_ref->{mbulge_count}=$1;
		}
		elsif($_=~m/mir-bulge range\t(\d+)__(\d+)/)
		{
			( $param_ref->{mbulge_begin},$param_ref->{mbulge_end} ) = ($1,$2);
		}
		elsif($_=~m/mir-bulge length\t(\d+)/)
		{
			$param_ref->{mbulge_len}=$1;
		}
		elsif($_=~m/utr-bulge allowed\t(\d+)/)
		{
			$param_ref->{ubulge_count} = $1;
		}
		elsif($_=~m/utr-bulge range\t(\d+)__(\d+)/)
		{
			($param_ref->{ubulge_begin},$param_ref->{ubulge_end} ) = ($1,$2);
		}
		elsif($_=~m/utr-bulge length\t(\d+)/)
		{
			$param_ref->{ubulge_len}=$1;
		}
		elsif($_=~m/loop allowed\t(\d+)/)
		{
			$param_ref->{loop_count} = $1;
		}
		elsif($_=~m/loop range\t(\d+)__(\d+)/)
		{
			( $param_ref->{loop_begin},$param_ref->{loop_end} ) = ($1,$2);
		}
		elsif($_=~m/loop length\t(\d+)/)
		{
			$param_ref->{loop_len}=$1;
		}
		elsif($_=~m/terminal mismatch\t(\d+)/)
		{
			$param_ref->{term_mismatch}=$1;
		}
		elsif($_=~m/loopORbulge allowed\t(\d+)/)
		{
			$param_ref->{mismatch_count}=$1;
		}
		elsif($_=~m/everything\t(\d+)/)
		{
			$param_ref->{all_count}=$1;
		}
		elsif($_=~m/loopORbulge range\t(\d+)__(\d+)/)
		{
			($param_ref->{mismatch_begin},$param_ref->{mismatch_end} ) = ($1,$2);
		}
	}
        
	$param_ref->{all_begin}  = $param_ref->{mismatch_begin};
	
	$param_ref->{all_end}    = $param_ref->{mismatch_end};
        
	close PARAM;
} 
	
	

sub scan_alignment
{
	my ( $seg1_utr_35, $seg2_utr_35, $seg3_mir_53, $seg4_mir_53, $feature_ref,$in,$source) = @_;
	my @seg1_utr_35 = split( //, $seg1_utr_35 );
	my @seg2_utr_35 = split( //, $seg2_utr_35 );
	my @seg3_mir_53 = split( //, $seg3_mir_53 );
	my @seg4_mir_53 = split( //, $seg4_mir_53 );
	my $matchcount = 0;
	# calculate term m,ismatch
	$feature_ref->{term_mismatch}       = cal_term_mismatch( $seg3_mir_53, $seg4_mir_53 );
	$feature_ref->{term_mismatch_mir35} = cal_term_mismatch_mir35( $seg3_mir_53, $seg4_mir_53 );
	$feature_ref->{matchseglen}         = cal_matchseglen($seg2_utr_35,$seg3_mir_53);
	my $align_len = length( $seg1_utr_35 );
	
	my $mir_pos = -1;
	
	my $utr_pos = -1;
	
	my $stem = 0; # only start counting features after the first stem
	
	my $loop_utropen = -1;
	
	#my $loop_utrclose = -1;
	
	my $loop_miropen = -1;
	
	#my $loop_mirclose = -1;
	
	my $ubulge_utropen = -1;
	
	#my $ubulge_utrclose = -1;
	
	my $ubulge_miropen = -1;
	
	#my $ubulge_mirclose = -1;
	
	my $mbulge_utropen = -1;
	
	#my $mbulge_utrclose = -1;
	
	my $mbulge_miropen = -1;
	
	#my $mbulge_mirclose = -1;
	
	
	for( my $i=0; $i<$align_len; $i++ ) # alignment pos is 0-based
	{
		# keep track of mir position
		if( $seg3_mir_53[$i] ne ' ' || $seg4_mir_53[$i] ne ' ' )
		{
			$mir_pos ++;
		}
		
		# keep track of utr position
		if( $seg1_utr_35[$i] ne ' ' || $seg2_utr_35[$i] ne ' ' )
		{
			$utr_pos ++;
		}
		

		
		##########################################################################################################
		# Going to insert chunk of codes for pending params below -for time being in $posBegin & $posEnd
		##########################################################################################################


		
		
		
		#GU pair
		if( ( $seg2_utr_35[$i] eq 'U' && $seg3_mir_53[$i] eq 'G' ) || ( $seg2_utr_35[$i] eq 'G' && $seg3_mir_53[$i] eq 'U' ) ||
				( $seg2_utr_35[$i] eq 'T' && $seg3_mir_53[$i] eq 'G' ) || ( $seg2_utr_35[$i] eq 'G' && $seg3_mir_53[$i] eq 'T' ) )
		{
			push @{$feature_ref->{gu}}, {	'utrpos' => $utr_pos,
																		'mirpos' => $mir_pos
																	};
		}				
		
		# set stem flag to 1 if match is found
		if($stem == 0)
		{
			if( $seg2_utr_35[$i] ne ' ' && $seg3_mir_53[$i] ne ' ' )
			{
				$stem = 1;	
			}
		}
		else
		{
			#stem
			
			if( $seg2_utr_35[$i] ne ' ' && $seg3_mir_53[$i] ne ' ' )
			{
				if( $loop_utropen == -1 && $ubulge_utropen == -1 && $mbulge_utropen == -1 )
				 {
				 }
				elsif( $loop_utropen > -1 )
				 {
					push @{$feature_ref->{loop}}, {	'utropen'  => $loop_utropen,
																					'utrclose' => $utr_pos,
																					'miropen'  => $loop_miropen,
																					'mirclose' => $mir_pos
																				};
																					
					$loop_utropen = - 1;
					
					$loop_miropen = - 1;
				}
				elsif( $ubulge_utropen > -1 )
				{
					push @{$feature_ref->{ubulge}}, {	'utropen'  => $ubulge_utropen,
																						'utrclose' => $utr_pos,
																						'miropen'  => $ubulge_miropen,
																						'mirclose' => $mir_pos
																					};
					
					$ubulge_utropen = -1;
					
					$ubulge_miropen = -1;
				}
				elsif( $mbulge_utropen > -1 )
				{
					push @{$feature_ref->{mbulge}}, {	'utropen'  => $mbulge_utropen,
																						'utrclose' => $utr_pos,
																						'miropen'  => $mbulge_miropen,
																						'mirclose' => $mir_pos
																					};
					
					$mbulge_utropen = -1;
					
					$mbulge_miropen = -1;
				}
			}
			# loop
			elsif( $seg1_utr_35[$i] ne ' ' && $seg4_mir_53[$i] ne ' ' )
			{
				if( $loop_utropen == -1 )
				{
					if(	$ubulge_utropen == -1 && $mbulge_utropen == -1 )
					{
						$loop_utropen = $utr_pos - 1;
						
						$loop_miropen = $mir_pos - 1;
					}
					elsif( $ubulge_utropen > -1)
					{
						$loop_utropen = $ubulge_utropen;
						
						$loop_miropen = $ubulge_miropen;
						
						$ubulge_utropen = -1;
						
						$ubulge_miropen = -1;
					}
					elsif( $mbulge_utropen > -1 )
					{
						$loop_utropen = $mbulge_utropen;
						
						$loop_miropen = $mbulge_miropen;
						
						$mbulge_utropen = -1;
						
						$mbulge_miropen = -1;
					}	
					
				}
				else
				{
				}
			}
			# mir bulge
			elsif( $seg1_utr_35[$i] eq ' ' && $seg2_utr_35[$i] eq ' ' )
			{
				if( $mbulge_utropen == -1 )
				{
					if( $loop_utropen > -1 )
					{
					}	
					else
					{	
						$mbulge_utropen = $utr_pos;
						
						$mbulge_miropen = $mir_pos - 1;
					}
				}
				else
				{
				}
			}
			# utr bulge	
			elsif( $seg3_mir_53[$i] eq ' ' && $seg4_mir_53[$i] eq ' ' )
			{
				if( $ubulge_utropen == -1)
				{
					if( $loop_utropen > -1 )
					{
					}	
					else
					{	
						$ubulge_utropen = $utr_pos - 1;
						
						$ubulge_miropen = $mir_pos;
					}
				}
				else
				{
				}
			}
		}
	}
	
	$feature_ref->{matchcount} = $matchcount;
}		
	

sub cal_term_mismatch
{
	my ( $seg3_mir_53, $seg4_mir_53 ) = @_;
	
	my ($blank_seg3) = $seg3_mir_53 =~ /^( *)/;
	
	my ($blank_seg4) = $seg4_mir_53 =~ /^( *)/;
	
	my $term_mismatch = length($blank_seg3) - length($blank_seg4);
	
	$term_mismatch    = 0 if $term_mismatch <0;
	#print"5'-Terminal mismatch:$term_mismatch\t";
	return $term_mismatch;
}


sub cal_term_mismatch_mir35
{
	my ( $seg3_mir_53, $seg4_mir_53 ) = @_;

	my $seg3_mir35   = reverse($seg3_mir_53);
	
	my $seg4_mir35   = reverse($seg4_mir_53);
	
	my ($blank_seg3) = $seg3_mir35 =~ /^( *)/;
	
	my ($blank_seg4) = $seg4_mir35 =~ /^( *)/;
	
	my $term_mismatch_mir35 = length($blank_seg3) - length($blank_seg4);
	
	$term_mismatch_mir35 = 0 if $term_mismatch_mir35 <0;
	return $term_mismatch_mir35;
}



sub cal_matchseglen
{
	my ( $seg2_utr_35, $seg3_mir_53 ) = @_;

	my $matchseglen = 0;

	my @seg2_utr_35 = split(//,$seg2_utr_35);

        my @seg3_mir_53 = split(//,$seg3_mir_53); 

	for( my $i=0; $i<scalar(@seg2_utr_35); $i++ )
	 {
         if($seg2_utr_35[$i] ne ' ' && $seg3_mir_53[$i] ne ' ')
		 {
                    $matchseglen++;
		 }
		 elsif($matchseglen>1 && $seg2_utr_35[$i] eq ' ')
		 {
			         last;
		 }
		 else
		 {
		 }
	 }
    return $matchseglen;
}





sub scan_alignment4misc_params
{
   my ($seg1_utr_35,$seg2_utr_35,$seg3_mir_53,$seg4_mir_53,$feature_ref,$in,$source) = @_;
   my $mir_pos=0;
   my $utr_pos=0;
   my @positions =('1_24','1_12','13_24','1_8','9_16','17_24');
   
   my @seg1_utr_35      = split(//,$seg1_utr_35);
   my @seg2_utr_35      = split(//,$seg2_utr_35);
   my @seg3_mir_53      = split(//,$seg3_mir_53);
   my @seg4_mir_53      = split(//,$seg4_mir_53);
   
   
   my $term_mismatch_3p = $feature_ref->{term_mismatch_mir35};
   my $term_mismatch_5p = $feature_ref->{term_mismatch};
      
      
   my $mir_len          = get_mirlength($seg3_mir_53,$seg4_mir_53); # Sub included
   
   my %feature_misc;
   for(my $j=0;$j<scalar(@positions);$j++)
		{
		   my ($posBegin,$posEnd)   = split(/_/,$positions[$j]);
		   my $matchcount           = 0;
		   my $mismatchbasecount    = 0;
		   my $mismatchutrcount     = 0;
		   my $gap                  = 0;
           
		   $mir_pos                 = 0;
		   $utr_pos                 = 0;
           for( my $i=0; $i<scalar(@seg1_utr_35); $i++ ) # alignment pos is 0-based
	        {
		       
		       if( $seg3_mir_53[$i] ne ' ' || $seg4_mir_53[$i] ne ' ' ) # keep track of mir position
		       {
			      $mir_pos ++;
		       }
				       
		       if( $seg1_utr_35[$i] ne ' ' || $seg2_utr_35[$i] ne ' ' ) # keep track of utr position
		       {
			      $utr_pos ++;
		       }

			   
		       if($mir_pos>=$posBegin && $mir_pos<=$posEnd&& $mir_pos<=($mir_len-$term_mismatch_3p) && $mir_pos>=$term_mismatch_5p)  # (1) Keep track of total matches in the site
		       {
		          if($seg3_mir_53[$i] ne ' ')
		          {
		              $matchcount++;
		          }
		       }
		
		                 
		      if($mir_pos>=$posBegin && $mir_pos<=$posEnd&&$mir_pos<=($mir_len-$term_mismatch_3p)&& $mir_pos>=$term_mismatch_5p) # (2) Keep track of total mismatches(bases wrt miRNA)
		      {
		         if($seg4_mir_53[$i] ne ' ')
		          {
		              $mismatchbasecount++;
		          }
 		      }
			
		                           
		     if($mir_pos>=$posBegin && $mir_pos<=$posEnd &&$mir_pos<=($mir_len-$term_mismatch_3p)&& $mir_pos>=$term_mismatch_5p)  # (3) Keep track of mismatches in UTR   
		     {
		       if($seg2_utr_35[$i] eq ' ' && $seg1_utr_35[$i] ne ' ')
		        {
		              $mismatchutrcount++;
		        }
		     }
			

			  if($mir_pos>=$posBegin && $mir_pos<=$posEnd &&$mir_pos<=($mir_len-$term_mismatch_3p)&& $mir_pos>=($term_mismatch_5p))  # (4) Keep track of internal gaps
			  {
			    if( $seg4_mir_53[$i] ne ' ' && $seg1_utr_35[$i] eq ' ')
			    {
				    $gap++;
			    }
		      }
		 }
        
		$feature_misc{'mir_length'}            = $mir_len;
		my $temp_matchcount                    = 'matchcount'.$posBegin.'_'.$posEnd;
        my $temp_mismatchbasecount             = 'mismatchbasecount'.$posBegin.'_'.$posEnd;
        my $temp_mismatchutrcount              = 'mismatchutrcount'.$posBegin.'_'.$posEnd;
		my $temp_gap                           = 'gap'.$posBegin.'_'.$posEnd;		
        $feature_misc{$temp_matchcount}        = $matchcount;
        $feature_misc{$temp_mismatchbasecount} = $mismatchbasecount;
        $feature_misc{$temp_mismatchutrcount}  = $mismatchutrcount; 
		$feature_misc{$temp_gap}               = $gap; 
   }
		return \%feature_misc;
}
    
	
	
sub get_mirlength
{
	my ($seg3_mir_53, $seg4_mir_53 ) = @_;
	
	my $mir_len;
	
	for (my $i=0; $i<length($seg3_mir_53); $i++)
	{
		if ( substr( $seg3_mir_53, $i, 1 ) ne ' ' || substr( $seg4_mir_53, $i, 1 ) ne ' ' )
		{
			$mir_len++;
		}
	}
	return $mir_len;
}	
	
	

 
