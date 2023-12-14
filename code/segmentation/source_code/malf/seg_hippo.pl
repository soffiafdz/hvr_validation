#!/usr/bin/env perl
#
##########################################################################
# do fusion-based segmentation:
#
# --- compute non-linear registrations/segmentations for each target
# for each side in {left, right}
#    for each of the 'n' best targets stored in 06-xcr
#       compute the non-linear registation between the target and subject  in 05-mh{L,R} space
#       transform the labels from the target onto the subject through the transformation
#       clean the labels to remove CSF voxels
#
# --- fuse all the target labels together.
# for each side in {left, right}
#   fuse 'n' label sets together to create the final label set.
#
# see inputs/outputs below
#
# Adapted by S.Fernandez (Aug-2022) from:
# Louis Collins, October 29, 2008
# McConnell Brain Imaging Centre,
# Montreal Neurological Institute,
# McGill University
##########################################################################

use strict;
use warnings "all";
use Getopt::Tabular;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use File::Glob ':glob';


####################################################################
my($Help, $Usage, $me);
my(@opt_table, %opt, @args);

$me = &basename($0);
%opt = (
	'verbose'   => 1,
	'clobber'   => 0,
	'keeptmp'   => 1,        # if 1, delete tmpdir by default
	'fake'      => 0,        # if 1, don't run commands; just print them.
	'startlevel' => 1,
	'stoplevel'  => 11,
	'tmpdir'     => "",
	'xfmdir'     => "",
	'labeldir'   => "",
	'modelSEG'   => "/ipl/ipl27/sfernandez/hvr_validation/libraries/malf_dorothee/models_seg",
	'modelDIR'   => "/ipl/ipl27/sfernandez/hvr_validation/libraries/malf_dorothee/models",
	'modelMRI'   => "/ipl/ipl27/sfernandez/hvr_validation/libraries/malf_dorothee/models_mri",
	'modelLAB'   => "/ipl/ipl27/sfernandez/hvr_validation/libraries/malf_dorothee/labels/reduced",
	'modelXFM'   => "/ipl/ipl27/sfernandez/hvr_validation/libraries/malf_dorothee/models_xfms",
);

$Help = <<HELP;
| $me does fusion-based segmentation of HC
| Problems or comments should be sent to: louis.collins\@mcgill.ca
HELP

$Usage = "Usage: $me [options] input_mri native_mask subj_id  tal.xfm  basename labels.mnc\n" .
"where:\n" .
"  inputs:\n".
"    input_mri   - native N3 corrected MRI\n".
"    native_mask - voxel-registered brain mask in native space\n".
"    subj_id     - string to identify subject (to avoid using as target)\n".
"    tal.xfm     - the native to stereotaxic space transformation (ie, where you want the labels)\n".
"   \n".
" outputs\n".
"    basename    - used for temp xfms and labels (example: ADNI.s1234.m12)\n".
"    labels.mnc  - output fused label set\n\n";


@opt_table = (
	["-verbose", "boolean", 0, \$opt{verbose}, "be verbose" ],
	["-clobber", "boolean", 0, \$opt{clobber}, "clobber existing files" ],
	["-keeptmp", "boolean", 0, \$opt{keeptmp},
		"keep data in temporary directory" ],
	["-fake", "boolean", 0, \$opt{fake}, "do a dry run, (echo cmds only)" ],
	["-stoplevel", "integer", 1, \$opt{stoplevel}, "stop at N fusions",
		'<int>'],
	["-startlevel", "integer", 1, \$opt{startlevel}, "start at M fusion",
		'<int>'],
	['-tmpdir', "string", 1, \$opt{tmpdir},
		"specify working directory [default defined in /tmp ]",
		"<string>"],
	['-xfmdir', "string", 1, \$opt{xfmdir},
		"specify tmp dir to store xfms+defs [default defined in /tmp ]",
		"<string>"],
	['-labeldir', "string", 1, \$opt{labeldir},
		"specify tmp dir to store labels [default defined in /tmp ]",
		"<string>"],
	['-modelDIR', "string", 1, \$opt{modelDIR},
		"specify alternate directory for stx MRIs [default:$opt{modelDIR}]",
		"<string>"],
	['-modelMRI', "string", 1, \$opt{modelMRI},
		"specify alternate directory for template MRIs (.5mm) [default:$opt{modelMRI}]",
		"<string>"],
	['-modelSEG', "string", 1, \$opt{modelSEG},
		"specify alternate directory for template  MRIs (1mm) [default:$opt{modelSEG}]",
		"<string>"],
	['-modelLAB', "string", 1, \$opt{modelLAB},
		"specify alternate directory for template labels [default:$opt{modelLAB}]",
		"<string>"],
	['-modelXFM', "string", 1, \$opt{modelXFM},
		"specify alternate directory for template transforms [default:$opt{modelXFM}]",
		"<string>"]
);


# Check arguments
&Getopt::Tabular::SetHelp($Help, $Usage);
&GetOptions (\@opt_table, \@ARGV) || exit 1;


die $Usage if(! ($#ARGV == 5));

# Get the arguments.

my $input_mri   = shift;    # inputs
my $natmask     = shift;
my $subj_id     = shift;
my $stx_xfm   = shift;
my $bname     = shift;      # outputs
my $labels    = shift;

print "input_mri   $input_mri     \n" if ($opt{verbose});
print "subj_id     $subj_id       \n" if ($opt{verbose});
print "stx_xfm     $stx_xfm       \n" if ($opt{verbose});
print "modelmridir $opt{modelMRI} \n" if ($opt{verbose});
print "modellabdir $opt{modelLAB} \n" if ($opt{verbose});
print "modelxfmdir $opt{modelXFM} \n" if ($opt{verbose});
print "bname       $bname         \n" if ($opt{verbose});
print "labels      $labels        \n" if ($opt{verbose});


# check files/ directories

die "$me: Couldn't find input file: $input_mri\n\n" if (!-e $input_mri);
die "$me: Couldn't find input file: $stx_xfm\n\n" if (!-e $stx_xfm);

die "$me: Couldn't find model xfm directory:   $opt{modelXFM}\n\n" if (!-e $opt{modelXFM});
die "$me: Couldn't find model label directory: $opt{modelLAB}\n\n" if (!-e $opt{modelLAB});
die "$me: Couldn't find model mri directory:   $opt{modelMRI}\n\n" if (!-e $opt{modelMRI});

die "$me: Couldn't find output xfm directory: $opt{xfmdir}\n\n"   if (($opt{xfmdir} ne "") && (!-e $opt{xfmdir}));
die "$me: Couldn't find output label directory: $opt{labeldir}\n\n" if (($opt{labeldir} ne "") && (!-e $opt{labeldir}));

die "$me: $labels exists, use -clobber to overwrite\n\n" if (-e $labels && !$opt{clobber});

# set constants

my $modeldir=$opt{modelDIR};
my $lowres_template="$modeldir/HC_icbm_avg_152_t1_tal_nlin_symmetric_VI.mnc";
my $hires_template="$modeldir/HC_icbm_avg_152_t1_tal_nlin_symmetric_VI-0.5mm.mnc";


my $model_leftxcorrmask="$modeldir/HC_icbm_avg_152_t1_tal_nlin_symmetric_VI_labels_left_xcorrmask.mnc";
my $model_rightxcorrmask="$modeldir/HC_icbm_avg_152_t1_tal_nlin_symmetric_VI_labels_right_xcorrmask.mnc";


# make tmpdir & define temporary files
my $tmpdir;
if ( $opt{tmpdir} eq "" ) {
	$tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => $opt{keeptmp} );
}
else {
	$tmpdir = $opt{tmpdir};
}

if ( $opt{xfmdir} eq "" ) { $opt{xfmdir} = "$tmpdir/${bname}_xfms" };
if ( $opt{labeldir} eq "" ) { $opt{labeldir} = "$tmpdir/${bname}_labels"};

mkdir $opt{xfmdir}   if (! -e $opt{xfmdir});
mkdir $opt{labeldir} if (! -e $opt{labeldir});

print "xfmdir      $opt{xfmdir}   \n" if ($opt{verbose});
print "labeldir    $opt{labeldir} \n" if ($opt{verbose});

# declare the file names needed to store data as we run
my $tmp_mhc_xfm     = "$tmpdir/${bname}_mhc.xfm"; # refined transforamtion from native  input data to stx template space, but only for the central region of the brain as defined by the average template (not the indiv templates).
my $tmp_mhL_xfm     = "$tmpdir/${bname}_mhL.xfm"; # further refined transformation, for the left medical temporal lobe (using the model_leftxcorrmask  defined above)
my $tmp_mhc_mhL_xfm = "$tmpdir/${bname}_mhc-concat-mhL.xfm"; # from native to local medial temp lobe for the left side.
my $tmp_mhR_xfm     = "$tmpdir/${bname}_mhR.xfm";
my $tmp_mhc_mhR_xfm = "$tmpdir/${bname}_mhc-concat-mhR.xfm";
my $tmp_mhc_mnc     = "$tmpdir/${bname}_mhc.mnc";
my $tmp_mhL_mnc     = "$tmpdir/${bname}_mhL.mnc";
my $tmp_mhR_mnc     = "$tmpdir/${bname}_mhR.mnc";
my $outfile_xcorrs  = "$tmpdir/${bname}_xcorrs.txt"; # will contain normalized mutual information vaues (no longer cross-correlation) values for each MTL zone, comparing input MRI to each of the template mnc volumes. (We will sort this to find the best n templates)

my $changed = 0;

# do local linear registration for the HC/MTL region

if ( ! (-e $tmp_mhc_xfm && -e $tmp_mhc_mnc) ||  $opt{clobber}) {
	&do_cmd("./code/malf/local_hc_reg","-clob", "-init_xfm", $stx_xfm, $input_mri,
		"-source_mask", $natmask, $tmp_mhc_xfm);
	&do_cmd("mincresample","-sinc","-clobber","-like", $lowres_template,
		"-transformation", $tmp_mhc_xfm, $input_mri, $tmp_mhc_mnc);
	$changed = 1;
}
else {
	print "$tmp_mhc_xfm and tmp_mhc_mnc exist already\n";
}

# do local linear registration to the model for the left HC region only.

if ( ! -e $tmp_mhL_xfm || ! -e $tmp_mhc_mhL_xfm || ! -e $tmp_mhL_mnc ||  $opt{clobber} ||  $changed ) {
	&do_cmd("minctracc","-step", 1, 1, 1, "-simplex", 2, "-lsq6", "-est_center",
		"-clobber","-tol", 0.0005, "-debug", "-source_mask",
		$model_leftxcorrmask, "-model_mask",  $model_leftxcorrmask,
		$tmp_mhc_mnc, $lowres_template, $tmp_mhL_xfm);

	if ( -e $tmp_mhc_mhL_xfm ) { unlink $tmp_mhc_mhL_xfm ; }

	&do_cmd("xfmconcat", $tmp_mhc_xfm, $tmp_mhL_xfm, $tmp_mhc_mhL_xfm);

	&do_cmd("mincresample","-sinc","-clobber","-like", $lowres_template,
		"-transformation", $tmp_mhc_mhL_xfm, $input_mri, $tmp_mhL_mnc);

	$changed = 1;
}
else {
	print "$tmp_mhL_xfm and $tmp_mhL_mnc exist already\n";
}

# do local linear registration to the model for the right HC region only.

if ( ! -e $tmp_mhR_xfm || ! -e $tmp_mhc_mhR_xfm || ! -e $tmp_mhR_mnc ||  $opt{clobber} || $changed ) {
	&do_cmd("minctracc","-step", 1, 1, 1, "-simplex", 2, "-lsq6",
		"-est_center", "-clobber","-tol", 0.0005, "-debug", "-source_mask",
		$model_rightxcorrmask, "-model_mask",  $model_rightxcorrmask,
		$tmp_mhc_mnc, $lowres_template, $tmp_mhR_xfm);

	if ( -e $tmp_mhc_mhR_xfm ) { unlink $tmp_mhc_mhL_xfm ; }
	&do_cmd("xfmconcat", $tmp_mhc_xfm, $tmp_mhR_xfm, $tmp_mhc_mhR_xfm);

	&do_cmd("mincresample","-sinc","-clobber","-like", $lowres_template,
		"-transformation", $tmp_mhc_mhR_xfm, $input_mri, $tmp_mhR_mnc);

	$changed = 1;
}
else {
	print "$tmp_mhR_xfm and $tmp_mhR_mnc exist already\n";
}

# build_xcorr_list for each of the models for the left and right sides

my $side;
my $m;
my $mutual_info;
my @models = split(/\n/, `ls $opt{modelLAB} | cut -d _ -f 2`);
print "Read " . @models . " model filenames\n" if $opt{verbose};

if ( -e $outfile_xcorrs && ! $opt{clobber} && ! $changed ) {
	print "$outfile_xcorrs exists, skipping to non-linear registrations\n";
}
else {
	open( OUTFILE_XCORRS, ">$outfile_xcorrs") or die "Cannot open output file: $!";
	foreach $side ("left","right") {
		my $this_template;
		my $mri;
		my @output;
		my @text;

		if ($side eq "left") {
			$this_template = "mhL"; # this label indicates the MTL left 'custom' registered data for this template volume for he left side.
			$mri = $tmp_mhL_mnc;
		}
		else {
			$this_template = "mhR";
			$mri = $tmp_mhR_mnc;
		}

		foreach $m (@models) {
			@output = `./code/malf/minctracc-w-nmi $opt{modelSEG}/${m}_$this_template.mnc $mri  -source_mask $opt{modelDIR}/HC_icbm_avg_152_t1_tal_nlin_symmetric_VI_labels_${side}_templatemask.mnc -identity -nmi -step 1 1 1 -blur_pdf 9 -simplex 1 -tol 0.01 -lsq6 -est_center -clob $tmpdir/tmp.xfm`;
			@text = grep('Final', @output);
			$mutual_info = $text[$#text];
			chomp $mutual_info;
			$mutual_info =~ s/Final objective function value = //;
			print                "$mutual_info $m.$this_template $side\n" if $opt{verbose};
			print OUTFILE_XCORRS "$mutual_info $m.$this_template $side\n";
		}
	}
	close (OUTFILE_XCORRS);
}


##############################################################################
#
#  OK... now we are ready to do the segmentations
#     for left and right sides
#        for the N best matched models (as defined in a sorted $outfile_xcorrs)
#           do non-linear reg between subject and model (within HC+AG mask)
#        fuse N labels
#     done


# blur the source w/1mm kernel

my $tmpN3blurred="$tmpdir/${bname}_source_n3_1mm";

if ( ! -e "${tmpN3blurred}_blur.mnc" || $opt{clobber} ) {
	&do_cmd( 'mincblur', '-clobber', '-fwhm', 1, $input_mri, $tmpN3blurred);
}

$tmpN3blurred = "${tmpN3blurred}_blur.mnc";

my $stx_source = "$tmpdir/${bname}_source_resampled_stx.mnc";

if ( ! -e $stx_source || $opt{clobber} ) {
	&do_cmd('mincresample', '-clob','-like', $lowres_template, $tmpN3blurred,
		$stx_source, '-transform', $stx_xfm);
}

# Reshape left xcorrmask from 0-255 to 0-1
my $reshaped_l_xcorrmask = "$tmpdir/$bname.reshaped_left_xcorrmask.mnc";
if ( ! -e $reshaped_l_xcorrmask ) {
	my $left_xcorrmask = "$opt{modelDIR}/HC_icbm_avg_152_t1_tal_nlin_symmetric_VI_labels_left_xcorrmask.mnc";
	&do_cmd('mincreshape', '-clobber', '-image_range', 0, 1,
		$left_xcorrmask, $reshaped_l_xcorrmask);
}

my @numbers = ();

for $side ("left", "right") {
	my $xfm = $tmp_mhc_mhL_xfm;
	$xfm = $tmp_mhc_mhR_xfm if ($side eq "right" );

	my $source = "$tmpdir/source_resampled_$side.mnc";
	&do_cmd('mincresample', '-clobber', '-like', $hires_template, $tmpN3blurred,
		$source, '-transform', $xfm);

	my $tailamount = $opt{stoplevel} - $opt{startlevel}+ 1;

	# get the 'n' best models, where best=most similar after local reg for that side
	my @models=`grep $side $outfile_xcorrs | grep -v $subj_id | sort -n | head -$opt{stoplevel}| tail -$tailamount `;
	my $count = 0;
	for $m (@models) {
		$count++;
		my $num = sprintf("%02d",$count); # note the upper limit of 99 labels to be fused...
		my ($score, $model) = split (" ",$m);
		(my $model_ = $model) =~ s/\./_/g;
		print "model = $model with score = $score\n";

		my $target = "$opt{modelMRI}/${model_}_1mm_blur.mnc";
		my $lin_xfm = "$tmpdir/$bname.$num.lin.$model.$side.xfm"; # (specific to side) lin from stx space input data  to template
		my $nlin1_xfm = "$tmpdir/$bname.$num.nlin1.$model.$side.xfm"; # low res non-lin def (look at -step below)
		my $nlin2_xfm = "$tmpdir/$bname.$num.nlin2.$model.$side.xfm"; # high res non-lin def

		my $this_result_xfm = "$opt{xfmdir}/$bname.$num.$model.$side.fusion.xfm";
		if (-e $this_result_xfm) {
			print "$this_result_xfm exists already, skipping estimation\n";
		}
		else {
### do I want to change the tolerance for linear reg? ie, smaller at 0.005 or 0.001?
			&do_cmd('./code/malf/minctracc-w-nmi', '-source_mask',
				"$opt{modelDIR}/HC_icbm_avg_152_t1_tal_nlin_symmetric_VI_labels_${side}_xcorrmask0.5mm.mnc",
				'-identity','-nmi','-blur_pdf', 9, '-simplex', 1, '-tol', 0.01,
				'-lsq6', '-step', 1, 1, 1, '-est_center', '-clob',
				$target, $source, $lin_xfm);

			&do_cmd('minctracc', '-nonlin', 'corrcoeff', '-tol', 0.0001,
				'-iterations', 15, '-max', 20, '-weight', 1, '-stiff', 1, '-simil', 0.3,
				'-step', 4, 4, 4, '-sub_lattice', 8, '-lattice_dia', 15, 15, 15,
				'-transformation', $lin_xfm, '-clobber', '-debug', '-source_mask',
				"$opt{modelDIR}/HC_icbm_avg_152_t1_tal_nlin_symmetric_VI_labels_${side}_xcorrmask0.5mm.mnc",
				$target, $source, $nlin1_xfm);

			$nlin2_xfm = $this_result_xfm;

			&do_cmd('minctracc', '-nonlin', 'corrcoeff', '-tol', 0.0001,
				'-iterations', 10, '-max', 20, '-weight', 1, '-stiff', 1, '-simil', 0.3,
				'-step', 2, 2, 2, '-sub_lattice', 7, '-lattice_dia', 10, 10, 10,
				'-transformation', $nlin1_xfm, '-clobber', '-debug', '-source_mask',
				"$opt{modelDIR}/HC_icbm_avg_152_t1_tal_nlin_symmetric_VI_labels_${side}_xcorrmask0.5mm.mnc",
				$target, $source, $nlin2_xfm);
		}

		# do seg for this xfm
		my $tmpxfm = "$tmpdir/$bname.$num.all.$model.$side.xfm";
		unlink $tmpxfm if (-e $tmpxfm);
		&do_cmd('xfmtool',"-invert:$opt{modelXFM}/${model_}-TO-mdl.xfm",
			$this_result_xfm, "-invert:$xfm", $stx_xfm, $tmpxfm);

		my $mod = `echo $model | cut -d . -f 1`;
		chomp $mod;

		my $tmpHC = "$tmpdir/$bname.$num.all.$model.$side.mnc";
		if (-e $tmpHC) {
			print "$tmpHC exists already\n";
		}
		else {
			#### check to see if the -like should not be the input MRI, because thats where we want the final labels.
			&do_cmd('mincresample', "$opt{modelLAB}/tal_${mod}_hp_csf_ld.mnc",
				'-nearest', '-clob', '-labels', '-like', $lowres_template,
				'-transformation', $tmpxfm, $tmpHC);
		}

		# clean up CSF voxels
		my $tmpHCclean = "$tmpdir/$bname.$num.allclean.$model.$side.mnc";
		if (-e $tmpHCclean) {
			print "$tmpHCclean exists already\n";
		}
		else {
			# remove voxels that are highly likely to be CSF.  Leave everything else as is...
			&do_cmd('./code/malf/pp_hc_xcorr_only',$tmpHC, $stx_source ,$tmpHCclean);
		}

		if ($side eq 'right') { # then both left and right are done for this model. put left + right together
			push @numbers, $num ;

			#my @left_labels  = <$tmpdir/$bname.$num.allclean.*.left.mnc>;
			#my @right_labels = <$tmpdir/$bname.$num.allclean.*.right.mnc>;

			# Don't use clean ones, keep CSF
			my @left_labels  = <$tmpdir/$bname.$num.all.*.left.mnc>;
			my @right_labels = <$tmpdir/$bname.$num.all.*.right.mnc>;

			die "too many left labels:".join(":",@left_labels)."\n" if (@left_labels > 1);
			die "too many right labels:".join(":",@right_labels)."\n" if (@left_labels > 1);

			# if we are on the left side (indicated by a nonzero left xcorrmask), then store left label.
			# otherise we are on the right, and need to store the right side label.
			# in both cases, store a 99 for background...
			&do_cmd('minccalc', '-clobber', '-discrete', '-expression',
				'A[0] > 0.5 ? (A[1]>0.5 ? A[1] : 99) : (A[2]>0.5 ? A[2] : 99)',
				$reshaped_l_xcorrmask, $left_labels[0], $right_labels[0],
				"$opt{labeldir}/$bname.$num.allclean.mnc");
		}
	} # for all models
} # for side

# now, time to fuse with voting technique.

my ($i, $num);
my $previous_sum;
my $labels_fused;
my $tmp_labels_fused;

my @labs = (99, 11, 12, 21, 22);
my @crisp_file_lab_args;

my @labels_to_fuse = ();

my $lab;

foreach $num (@numbers) {

	$i = 1.0 * $num;

	my $n = sprintf("%02d",$num);
	my $label_set = "$opt{labeldir}/$bname.$n.allclean.mnc";


	push @labels_to_fuse, $label_set;

# average the labels for the 1..$num set of labels
	my @crisp_file_lab_args = ();
	foreach $lab (@labs) { # extract the label and create an average
		my $avg = "$tmpdir/$bname.$num.average.$lab.mnc";
		&do_cmd('mincaverage', '-clobber', '-binarize', '-binvalue', $lab, @labels_to_fuse, $avg);

		push @crisp_file_lab_args, $lab;
		push @crisp_file_lab_args, $avg;

	}

  # make a 'crisp' set of labels, essentially use max vote at each voxel
  # note that there is a slight bias, when two labels are tied (say both at 50%) the higher number label will win.
	$tmp_labels_fused = "$tmpdir/$bname.$num.fusion.mnc";
	#&do_cmd('crispify', '-clobber', '-volume', $tmp_labels_fused, @crisp_file_lab_args);
	&do_cmd('itk_merge_labels', '--byte', '--clobber',
		@crisp_file_lab_args, $tmp_labels_fused);

  # get rid of the background value '99' used above in the segmentation.
	$labels_fused = "$opt{labeldir}/$bname.$num.fusion.mnc";
	&do_cmd('minclookup', '-clobber', '-discrete', '-int',
		'-lut', '11 11; 12 12; 21 21; 22 22',
		$tmp_labels_fused, $labels_fused);
}

# copy last labels fused as result
&do_cmd('minccalc','-clobber', '-discrete', '-expression', 'A[0]',
	$labels_fused, $labels);

exit 0;

sub do_cmd {
	print STDOUT "@_\n" if $opt{verbose};
	if(!$opt{fake}){
		system(@_) == 0 or die;
	}
}
