#!/usr/bin/perl

use warnings; use strict; use feature ':5.10';
use Bio::DB::Fasta;
use Carp;
use Cwd;
use Data::Dumper;
use File::Basename;
use File::Path qw(make_path);
use File::Slurp qw/read_file/;
use File::Temp;
use List::MoreUtils qw/natatime/;
use Parse::RecDescent;
use Try::Tiny;

croak "Provide source directory" unless ($ARGV[0] and -d $ARGV[0]);
croak "Provide destination directory" unless ($ARGV[1] and -d $ARGV[1]);

execute($ARGV[0], $ARGV[1], "GRCh37");

sub execute {
    my ($root, $destination, $major) = @_;
    
    print "\nThis script requires several gigabytes of free space in /tmp.\n";
    print "Also, please ensure that the destination directory $destination is empty before continuing.\n";
    print "Any files in this directory may be overwritten or modified.\n";
    print "Continue (yes/no): ";
    chomp(my $res = <STDIN>);
    return unless $res eq "yes";

    my $store;

    $store->{alt_loci}  = {};
    $store->{patches}   = {};
    $store->{primary}   = {};
    $store->{flank}     = ($ARGV[2] and $ARGV[2] =~ /^\d+$/) ? $ARGV[2] : 150;
    $store->{destdir}   = getcwd() . '/' . $destination;
    $store->{tempdir}   = File::Temp->newdir();

    say "Destination: " . $store->{destdir};
    say "Temporary fasta storage: " . $store->{tempdir};
    say "Flank is $store->{flank}.";

    make_path($destination);

    $store->{reference} = grab_lite_reference($store, $root);
    
    $root =~ s/(.+)\/$/$1/; # trim off trailing /
    chdir $root;
    $root = ".";
    
    opendir(my $dh, $root) || croak "Can't opendir $root: $!";
    my @minors = sort {$a cmp $b} grep { $_ =~ /^$major(?:\.p\d+)?$/ and -d "$root/$_" } readdir($dh);
    closedir $dh;
    
    for my $minor (@minors) {
        process_release($store, "$root/$minor", $major, $minor);
    }
    
    my $alt_loci = $store->{alt_loci};
    my $patches  = $store->{patches};

    # process alt ref loci
    for my $rel (sort {$a cmp $b} keys %$alt_loci) {
        for my $assembly (sort {$a cmp $b} keys %{$alt_loci->{$rel}}) {
            my $assemblyname = sprintf "%s_%s", $rel, $assembly;
            printf "Processing %s\n", $assemblyname;
            create_remap($store, $alt_loci->{$rel}{$assembly}, $assemblyname);
            # create_remap calls process_segs, which returns nothing but can throw exceptions
            # exceptions are caught in create_remap; create_remap has no need to return anything
        }
    }
    
    # process patches
    for my $rel (sort {$a cmp $b} keys %$patches) {
        my $assemblyname = sprintf "%s_%s", $rel, 'PATCHES';
        printf "Processing %s\n", $assemblyname;
        create_remap($store, $patches->{$rel}, $assemblyname);
        # create_remap calls process_segs, which returns nothing but can throw exceptions
        # exceptions are caught in create_remap; create_remap has no need to return anything
    }
    
    say "Done!";
}

# extracts the lite reference and pulls seqids from it 
sub grab_lite_reference {
    my ($store, $root) = @_;

    my $tempdir  = $store->{tempdir};

    my $unzipped = "$tempdir/GRCh37-lite.fa";

    if (-e "/tmp/GRCh37-lite.fa") {
        say "Using /tmp/GRCh37-lite.fa";
        $unzipped = "/tmp/GRCh37-lite.fa";
    } else {
        my $gzipped  = "$root/GRCh37/special_requests/GRCh37-lite.fa.gz";

        say "Extracting $gzipped to $unzipped.";
        croak "GRCh37 lite not in expected place '$gzipped'" unless -e $gzipped;

        say "Skip this step by extracting GRCh37/special_requests/GRCh37-lite.fa.gz to /tmp/GRCh37-lite.fa before running this script.";

        my $cmd = "zcat $gzipped > $unzipped";
        system "$cmd";
    }
    
    open my $fafh, '-|', "grep '>' $unzipped" || croak "Unable to seq ids from reference fasta.";
    
    my @seqs;
    
    while (my $line = <$fafh>) {
        chomp $line;
        $line =~ /^
                >                   # line should begin with >
                (\S+)               # the seqid
                \s+                 # space is a seperator
                (.+)                # the description
            $/x || croak "Can't parse seqid '$line' in fasta";

            my $seqid = $1;
            my $desc  = $2;
            my $acc;

            if ($seqid =~ /(\d+|[XY])/) {
                $desc =~ /^
                        (
                            [a-zA-Z0-9]+    # acc is something like CM000663
                            (?:\.\d+)?      # with an optional suffix of the form .11
                        )
                        \s.+
                    $/x || croak "Can't parse accession out of descritpion '$desc'";
                $acc = $1;
            } else {
                $acc = $seqid;
            }

        my $fasta_seq = {
            acc   => $acc,
            seqid => $seqid,
            line  => $line,
            fasta => $unzipped,
        };
        $fasta_seq->{desc} = $desc || '';
        push @seqs, $fasta_seq;
    }
    
    close $fafh;

    return {$unzipped => \@seqs};
}

# processes a release of GRCh37
sub process_release {
    my ($store, $root, $major, $minor) = @_;
    
    my $ispatch = 0;
    
    if ($minor eq $major) {
        say "Processing base version $minor";
    } elsif ($minor =~ /^$major\.p(\d+)$/) {
        $ispatch = 1;
        say "Processing version $major patch $1";
    }
    
    opendir(my $dh, "$root") || croak "Can't opendir $root: $!";
    my @dirs = grep {$_ ne '.' and $_ ne '..' and -d "$root/$_"} readdir($dh);
    closedir $dh;
    
    my @locis     = sort {$a cmp $b} grep {$_ =~ /^ALT_REF_LOCI_\d+$/} @dirs;
    my @primaries = grep {$_ eq "Primary_Assembly"} @dirs;
    my @patches   = grep {$_ eq "PATCHES"} @dirs;
    
    croak "Did not find 9 ALT_REF_LOCI dirs" unless scalar(@locis) == 9;
    croak "Did not find Primary_Assembly dir" unless scalar(@primaries) == 1;
    croak "Did not find PATCHES dir" unless (!$ispatch or scalar(@patches) == 1);
    
    for my $primary (@primaries) {
        if (-l "$root/$primary") {
            my $link = readlink("$root/$primary");
            if ($link =~ /\.\.\/($major(?:\.p\d+)?)\/Primary_Assembly/) {
                if (defined $store->{primary}{$1}) {
                    $store->{primary}{$minor} = $store->{primary}{$1};
                } else {
                    croak "Primary_Assembly pointed to release not yet processed, must have done something out of order.";
                }
            } else {
                croak "Can't understand what link $primary -> $link refers to";
            }
        } else {
            $store->{primary}{$minor} = process_primary_assembly($store, "$root/$primary");
        }
    }
    
    for my $loci (@locis) {
        $loci =~ /^ALT_REF_LOCI_(\d+)$/ || croak "Can't parse $loci";
        my $num = $1;
        
        if (-l "$root/$loci") {
            my $link = readlink("$root/$loci");
            if ($link =~ /\.\.\/($major(?:\.p\d+)?)\/(ALT_REF_LOCI_\d+)/) {
                croak "Loci $loci links to $2 but $loci != $2" if $2 ne $loci;
                if (defined $store->{alt_loci}{$1}{$loci}) {
                    $store->{alt_loci}{$minor}{$loci} = $store->{alt_loci}{$1}{$loci};
                } else {
                    croak "ALT LOCI pointed to release not yet processed, must have done something out of order.";
                }
            } else {
                croak "Can't understand what link $loci -> $link refers to";
            }
        } else {
            $store->{alt_loci}{$minor}{$loci} = process_alt_loci($store, "$root/$loci");
        }
        
    }
    
    for my $patch (@patches) {
        $store->{patches}{$minor} = process_patches($store, "$root/$patch");
    }
}

# creates a mapping between chromosomes and accessions in the primary reference
sub process_primary_assembly {
    my ($store, $root) = @_;
    
    my $chr2acc = "$root/assembled_chromosomes/chr2acc";
    my %acc2chr;
    
    croak "Can't read required file $chr2acc from the primary assembly" unless -e $chr2acc;
    
    open my $chr2accfh, '<', $chr2acc;
    
    while (my $line = <$chr2accfh>) {
        chomp $line;
        $line =~ /^
            (
                \d{1,2}|X|Y
            ) \t (
                [a-zA-Z0-9]+
                (?:\.\d+)
            )
        $/x || croak "Can't parse '$line' from $chr2acc";
        
        my $chr = $1;
        my $acc = $2;
        
        croak "Duplicate acc '$acc' found in $chr2acc" if exists $acc2chr{$acc};
        
        $acc2chr{$acc} = $chr;
    }
    
    close $chr2accfh;
    return {acc2chr => \%acc2chr};
}

# finds the scaffolds in an ALT_REF_LOCI_N directory
sub process_alt_loci {
    my ($store, $root) = @_;
    
    opendir(my $dh, "$root") || croak "Can't opendir $root: $!";
    my @dirs = sort {$a cmp $b} grep {$_ ne '.' and $_ ne '..' and -d "$root/$_" and "$root/$_" =~ /(?:alt|placed)_scaffolds/} readdir($dh);
    closedir $dh;
    
    croak "Did not find expected directory in $root" unless scalar(@dirs) == 1;
    
    return process_scaffolds($store, "$root/$dirs[0]");
}

# finds the scaffolds in a PATCHES directory
sub process_patches {
    my ($store, $root) = @_;
    
    opendir(my $dh, "$root") || croak "Can't opendir $root: $!";
    my @dirs = sort {$a cmp $b} grep {$_ ne '.' and $_ ne '..' and -d "$root/$_" and "$root/$_" =~ /(?:alt|placed)_scaffolds/} readdir($dh);
    closedir $dh;
    
    croak "Did not find expected directory in $root" unless scalar(@dirs) == 1;
    
    return process_scaffolds($store, "$root/$dirs[0]");
}

# proceses the scaffolds directory
sub process_scaffolds {
    my ($store, $root) = @_;
    
    croak "Can't find FASTA directory in $root" unless -d "$root/FASTA";
    croak "Can't find alignments directory in $root" unless -d "$root/alignments";
    
    my @gzfastas = @{gather_fasta($store, "$root/FASTA")};
    my @asns     = @{gather_asn($store, "$root/alignments")};
    
    croak "More than one FASTA file found in $root, not sure what to do" unless scalar(@gzfastas) == 1;
    
    # zcat all fasta files
    my $fastas     = process_fastas($store, $root, \@gzfastas);
    my $patches    = process_asns($store, $root, \@asns);
    my $placements = process_placements($store, $root);
    
    # join patches and placements together
    # TODO check this
    for (keys %$patches) {
        if (exists $placements->{$_}) {
            $patches->{$_}{placement} = $placements->{$_};
            delete $placements->{$_};
        } else {
            # try harder to match accessions by trimming off the suffix number
            # create a mapping of trimmed accs to untrimmed accs
            my %acc_remap;
            for (keys %$placements) {
                $_ =~ /(.+)\.\d+-(.+)\.\d+$/ || croak "Couldn't parse accession $_";
                croak "There are multiple versions of acc '$1-$2'" if (exists $acc_remap{"$1-$2"}); 
                $acc_remap{"$1-$2"} = $_;
            }
            
            if (exists $acc_remap{$_}) {
                my $full_acc = $acc_remap{$_};
                $patches->{$_}{placement} = $placements->{$full_acc};
                delete $placements->{$full_acc};
            } else {
                warn colorize('yellow', "Could not find placement for patch '$_'\n");
            }
        }
    }
    croak "Did not find patch for all placements" unless (scalar keys %$placements) == 0;
    
    return {
        fastas     => $fastas,
        alignments => $patches,
    };
}

# gathers fasta files in a directory
sub gather_fasta {
    my ($store, $root) = @_;
    
    opendir(my $dh, "$root") || croak "Can't opendir $root: $!";
    my @dirs = sort {$a cmp $b} grep {$_ ne '.' and $_ ne '..' and $_ =~ /.+\.fa.gz$/} readdir($dh);
    closedir $dh;
    
    return \@dirs;
}

# gathers asn files in a directory
sub gather_asn {
    my ($store, $root) = @_;
    
    opendir(my $dh, "$root") || croak "Can't opendir $root: $!";
    my @dirs = sort {$a cmp $b} grep {$_ ne '.' and $_ ne '..' and $_ =~ /.+\.asn$/} readdir($dh);
    closedir $dh;
    
    return \@dirs;
}

# extracts fastas and pulls seqids
sub process_fastas {
    my ($store, $root, $gzfastas) = @_;
    
    my $tempdir = $store->{tempdir};
    
    my %fastas;
    
    for my $fa (@$gzfastas) {
        $fa =~ /(.+fa)\.gz$/ || croak "Can't parse fasta file name $fa";
        my $newfa = $1;
        
        $root =~ /^\.\/(.+)$/ || croak "Can't trim leading . from $root";
        my $newroot = $1;
        
        my $gzipped      = "$root/FASTA/$fa";
        my $unzipped_dir = "$tempdir/$newroot/FASTA";
        my $unzipped     = "$unzipped_dir/$newfa";

        make_path("$unzipped_dir");
        
        my $cmd = "zcat $gzipped > $unzipped";
        # say "\t\t$cmd";
        system "$cmd";
        
        # pull out the accession and gi numbers in this fasta file
        open my $fafh, '-|', "grep '>' $unzipped" || croak "Unable to pull acc and gi nums from fasta file.";
        
        my @seqs;
        
        while (my $line = <$fafh>) {
            chomp $line;
            $line =~ /^
                    > gi \|             # line should begin with >gi|
                    (
                        \d+             # gi is a number
                    )
                    \| gb \|            # then, should encounter |gb|
                    (
                        [a-zA-Z0-9]+    # acc is something like GL949741
                        (?:\.\d+)?      # with an optional suffix of the form .11
                    )
                    \|                  # then a final | in the seqid
                    (\s.+)?             # optional whitespace and description after seqid
                $/x || croak "Can't parse seqid '$line' in fasta";
                
            my $fasta_seq = {
                gi    => $1,
                acc   => $2,
                seqid => "gi|$1|gb|$2|",
                line  => $line,
                fasta => $unzipped,
            };
            $fasta_seq->{desc} = $3 || '';
            push @seqs, $fasta_seq;
        }
        
        close $fafh;
        $fastas{$unzipped} = \@seqs;;
    }
    
    return \%fastas;
}

# processes asn files but does not parse data in them
sub process_asns {
    my ($store, $root, $asns) = @_;
    
    my %patches;
    
    for my $asn (@$asns) {
        $asn =~ /
                ^
                (
                    [a-zA-Z0-9]+    # acc is something like GL949741
                    (?:\.\d+)?      # with an optional suffix of the form .11
                )
                _
                (
                    [a-zA-Z0-9]+    # acc is something like GL949741
                    (?:\.\d+)?      # with an optional suffix of the form .11
                )
                \.asn$
            /x || croak "Can't parse asn file name '$asn'";
        
            
        my $src_acc = $1;
        my $dst_acc = $2;
        
        my $gff = sprintf("%s_%s.gff", $src_acc, $dst_acc);

        croak "ASN '$root/alignments/$asn' did not exist alongside gff" unless (-s "$root/alignments/$asn");
        croak "GFF '$root/alignments/$gff' did not exist alongside asn" unless (-s "$root/alignments/$gff");
        croak "Duplicate accession '$src_acc-$dst_acc' found in asn/gff dir" if exists $patches{"$src_acc-$dst_acc"};
        
        $patches{"$src_acc-$dst_acc"} = {
            src => "$src_acc",
            dst => "$dst_acc",
            asn => "$root/alignments/$asn",
            gff => "$root/alignments/$gff",
        };
    }
    
    return \%patches;
}

# parses placement file
sub process_placements {
    my ($store, $root) = @_;
    
    my %placements;
    
    open my $placementfh, '<', "$root/alt_scaffold_placement.txt" || croak "Couldn't open $root/alt_scaffold_placement.txt";
    
    while (my $line = <$placementfh>) {
        chomp $line;
        if ($line =~ /^#/) {
            croak "Unrecognized header in alt_scaffolds_placement.txt" if 
                $line ne "#alt_asm_name\tprim_asm_name\talt_scaf_name\talt_scaf_acc\t".
                "parent_type\tparent_name\tparent_acc\tregion_name\tori\talt_scaf_start\t".
                "alt_scaf_stop\tparent_start\tparent_stop\talt_start_tail\talt_stop_tail";
            next;
        }
        
        my @fields = split /\t/, $line;
        my $src_acc = $fields[3];
        my $dst_acc = $fields[6];
        
        croak "Duplicate accession '$src_acc-$dst_acc' found in alt_scaffold_placement.txt" if exists $placements{"$src_acc-$dst_acc"};

        $placements{"$src_acc-$dst_acc"} = {
            alt_asm_name   => $fields[0],
            prim_asm_name  => $fields[1],
            alt_scaf_name  => $fields[2],
            alt_scaf_acc   => $fields[3],
            parent_type    => $fields[4],
            parent_name    => $fields[5],
            parent_acc     => $fields[6],
            region_name    => $fields[7],
            ori            => $fields[8],
            alt_scaf_start => $fields[9],
            alt_scaf_stop  => $fields[10],
            parent_start   => $fields[11],
            parent_stop    => $fields[12],
            alt_start_tail => $fields[13],
            alt_stop_tail  => $fields[14],
        };
    }
    
    return \%placements;
}

# for verifying accessions
sub alignment_get_acc {
    my ($store, $alignment) = @_;
    
    if (exists $alignment->{placement}) {
        return (
            $alignment->{placement}{alt_scaf_acc},
            $alignment->{placement}{parent_acc}
        );
    }
    return (
        $alignment->{src},
        $alignment->{dst}
    );
}

# runs through patches and alternates and prepares for remap file creation
sub create_remap {
    my ($store, $scaffolds, $outputname) = @_;
    # outputname refers to the remap file to write to
    
    my $fastas     = $scaffolds->{fastas};
    my $alignments = $scaffolds->{alignments};
    
    for my $alignment (values %$alignments) {

        my ($srcacc, $dstacc) = alignment_get_acc($store, $alignment);
        
        unless (exists $alignment->{remap}) {
            try {
                # get a fasta sequence hash given an accession
                # the hash contains the source .fa the seqid appears in
                my $par_fasta_seq = find_fasta_from_acc($store, $store->{reference}, $dstacc);
                my $alt_fasta_seq = find_fasta_from_acc($store, $fastas, $srcacc);
            
                my $seqalign = parse_asn($store, $alignment->{asn});
                
                # for now write straight to a file and record that it didn't throw an exception
                process_segs( {
                    store         => $store,
                    alignment     => $alignment,
                    par_fasta_seq => $par_fasta_seq,
                    alt_fasta_seq => $alt_fasta_seq,
                    seqalign      => $seqalign,
                    outputname    => $outputname
                } );
                $alignment->{remap} = {
                    status => 'success',
                    name   => $outputname,
                    par_fa => $par_fasta_seq,
                    alt_fa => $alt_fasta_seq,
                };
                say colorize('blue', "Processed $alignment->{asn}");
                1;
            } catch {
                chomp $_;
                print colorize('red');
                say "Parse error: $_"; # not $@
                say "        asn: $alignment->{asn}\n";
                print colorize('clear');
                $alignment->{remap} = {
                    status => 'failure',
                    error  => $_,
                    #par_fa => $par_fasta_seq,
                    #alt_fa => $alt_fasta_seq,
                };
                print Dumper($alignment);
            };
        }
    }
}

# helps match a fasta and accession
sub find_fasta_from_acc {
    my ($store, $fastas, $acc) = @_;
    
    my @found_files;
    
    # for each fasta file
    for my $fa (keys %$fastas) {
        # grep through all accessions for this single fasta
        my @matches = grep { $_->{acc} eq $acc } @{$fastas->{$fa}};
        if (scalar(@matches) == 1) {
            push @found_files, $matches[0];
        } elsif (scalar(@matches) == 0) {
            next;
        } else {
            croak "Found accession '$acc' multiple times in '$fa'";
        }
    }
    
    if (scalar(@found_files) != 1) {
        if (scalar(@found_files) == 0) {
            croak "Did not find accession '$acc' in any fasta files";
        }
        croak "Found accession '$acc' in multiple fasta files";
    }
    return $found_files[0];
}

# processes segs from an asn file
sub process_segs {
    my $args = shift;

    my $store          = $args->{store};
    my $alignment     = $args->{alignment};
    my $par_fasta_seq = $args->{par_fasta_seq};
    my $alt_fasta_seq = $args->{alt_fasta_seq};
    my $seqalign      = $args->{seqalign};
    my $outputname    = $args->{outputname};

    my $asn       = $alignment->{asn};
    my $placement = $alignment->{placement} if defined $alignment->{placement};
    my $par_fasta = $par_fasta_seq->{fasta};
    my $par_seqid = $par_fasta_seq->{seqid};
    my $alt_fasta = $alt_fasta_seq->{fasta};
    my $alt_seqid = $alt_fasta_seq->{seqid};

    my $chromosome = parse_chromosome($store, $alignment, $placement);
    
    my $par_fa = new Bio::DB::Fasta($par_fasta);
    my $alt_fa = new Bio::DB::Fasta($alt_fasta);
    
    # Make sure the id is actually there.
    my @par_full_ids    = $par_fa->ids();
    my @alt_full_ids    = $alt_fa->ids();
    croak "Couldn't narrow down to a single seqid in the ref fasta given seqid '$par_seqid'"
        unless scalar(grep {$par_seqid eq $_} @par_full_ids) == 1;
    croak "Couldn't narrow down to a single seqid in the alt fasta given seqid '$alt_seqid'"
        unless scalar(grep {$alt_seqid eq $_} @alt_full_ids) == 1;

    my $outfasta  = sprintf "%s/%s.fa", $store->{destdir}, $outputname;
    my $outremap  = sprintf "%s/%s.remap", $store->{destdir}, $outputname;
    
    open my $outfastafh, '>>', $outfasta || croak "Couldn't open $outfasta for writing";
    open my $outremapfh, '>>', $outremap || croak "Couldn't open $outremap for writing";
    
    my $count = 0;
    my $segs = get_seg_list($store, $seqalign->{type}, $seqalign->{segs});

    for my $subsegs (@{$segs}) {
        # Get ready to process the parsed asn file
        my $ids       = $subsegs->{ids};
        my $starts    = $subsegs->{starts};
        my $lens      = $subsegs->{lens};
        my $strands   = $subsegs->{strands};
        my $dim       = $subsegs->{dim};
        my $numseg    = $subsegs->{numseg};
        my $revstrand = 0;
        
        my $ori = parse_orientation($store, $strands, $placement);

        my @alt;
        my @par;
        my @len;
        my @ops;
        my @running_length;

        # TODO This assumes that the alt scaffold is first in the pair and the parent is second in the pair
        # TODO Does this assumption have consequences for 'negative' strand handling? 
        my $aa = 0;
        my $pp = 1;

        for (my $i = 0; $i < $numseg; $i++) {
            my @spos   = @{$starts->[$i]};
            my $seglen = $lens->[$i];

            my $op;
            if (scalar(grep {$_ == -1} @spos) == 0) {
                $op = 'M';
                push @alt, {
                    start => $spos[$aa],
                    stop  => $spos[$aa] + $seglen,
                };
                push @par, {
                    start => $spos[$pp],
                    stop  => $spos[$pp] + $seglen,
                };
            } elsif ($spos[$aa] == -1) {
                $op = 'D';
                push @par, {
                    start => $spos[$pp],
                    stop  => $spos[$pp] + $seglen,
                };
            } elsif ($spos[$pp] == -1) {
                $op = 'I';
                push @alt, {
                    start => $spos[$aa],
                    stop  => $spos[$aa] + $seglen,
                };
            } else {
                croak "Nonsense starting positions at index $i: " . join(', ', @spos);
            }

            push @len, $seglen;
            push @ops, $op;

            if ($op ne 'D') {
                my $segseq = $alt_fa->seq($alt_seqid, $alt[$#alt]->{start}+1 => $alt[$#alt]->{stop});
                
                if ($segseq =~ /^(N+)$/) {
                    my $nlen = length($1);
                    croak "Expected an insertion during a split (since N in alt should not map to parent)" if $op ne 'I';
                    croak "Previous cigar op was not a match during a split" if $ops[$#ops-1] ne 'M';

                    # The purpose of this is to remove the current entry when encountering
                    # a split. Since this only occurs during an 'I' the current entry will not
                    # have been added to @par, and there is no need to pop from @par.
                    pop @alt;
                    pop @len;
                    pop @ops;

                    process_remap_chunk( {
                        store         => $store,
                        ori           => $ori,
                        alt           => \@alt,
                        par           => \@par,
                        len           => \@len,
                        ops           => \@ops,
                        par_fasta_seq => $par_fasta_seq,
                        alt_fasta_seq => $alt_fasta_seq,
                        par_fa        => $par_fa,
                        alt_fa        => $alt_fa,
                        count         => $count,
                        chromosome    => $chromosome,
                        outremapfh    => $outremapfh,
                        outfastafh    => $outfastafh,
                    } );

                    @alt = ();
                    @par = ();
                    @len = ();
                    @ops = ();

                    $count++; 
                    $i++; # look at the next seg
                    @spos = @{$starts->[$i]};
                    croak "Excised a sequence of length $nlen (I) but it was not followed by a D."
                        unless ($spos[0] == -1 and $op eq 'I');
                    my $difflen = $nlen - $lens->[$i];
                    #warn colorize('gray', sprintf( "Excised N sequence was of ".
                    #    "length %6d but following D was only of length %6d; diff: %d. This is probably harmless.\n",
                    #    $nlen, $lens->[$i], $difflen) unless ($difflen eq 0));
                    next;
                }
            }
        }
        process_remap_chunk( {
            store         => $store,
            ori           => $ori,
            alt           => \@alt,
            par           => \@par,
            len           => \@len,
            ops           => \@ops,
            par_fasta_seq => $par_fasta_seq,
            alt_fasta_seq => $alt_fasta_seq,
            par_fa        => $par_fa,
            alt_fa        => $alt_fa,
            count         => $count,
            chromosome    => $chromosome,
            outremapfh    => $outremapfh,
            outfastafh    => $outfastafh,
        } );
        $count++;
    }

    close $outfastafh;
    close $outremapfh;
}

# makes sense of what chromosome a patch/alternate is referring to on the primary reference
sub parse_chromosome {
    my ($store, $alignment, $placement) = @_;

    my $chromosome;
    my $acc2chr = $store->{primary}{GRCh37}{acc2chr};
    if (defined $placement) { # TODO clean this up
        $chromosome = $placement->{parent_name};
        
        # Make sure chromosome matches up.
        # TODO should use the current release instead of just assuming GRCh37
        croak "Chromosome name specified by placements does not match acc2chr in primary assembly." unless
            (exists $acc2chr->{$placement->{parent_acc}} and $acc2chr->{$placement->{parent_acc}} eq $chromosome);
    } else {
        $chromosome = $acc2chr->{$alignment->{dst}};
    }
    croak "Couldn't find chromosome name $alignment->{dst}, this usually means this isn't mapping back onto something in the primary assembly" unless defined $chromosome; 
    
    return $chromosome;
}

# attempts to make sense of the orientation from placements file and patch file
sub parse_orientation {
    my ($store, $strands, $placement) = @_;

    # Verify that strands are consistent
    my $one = $strands->[0][0];
    my $two = $strands->[0][1];
    for (@$strands) {
        if ($_->[0] ne $one and $_->[1] ne $two) {
            croak "Strands changed direction among same sequence; not sure what this means. (Expected $one and $two but got " . $_->[0] . " and " . $_->[1] . ".)";
        }
    }

    # Figure out what the strands and orientation are saying
    # TODO some of these are unsupported, as they've yet to be encountered
    my $ori;
    if (defined $placement) {
        my %recognized_ori = (
            plus => {
                plus => { '+' => 'positive', '-' => 'invalid', 'b' => 'positive', },
                minus => { '+' => 'invalid', '-' => 'unsupported', 'b' => 'unsupported', },
            },
            minus => {
                plus => { '+' => 'invalid', '-' => 'negative', 'b' => 'negative', },
                minus => { '+' => 'unsupported', '-' => 'invalid', 'b' => 'unsupported', },
            },
        );
        
        $ori = $recognized_ori{$one}{$two}{$placement->{ori}} || croak "Orientation $one - $two: $placement->{ori} is not recognized";
        croak "Orientation $one - $two: $placement->{ori} is $ori" unless grep {$ori eq $_} qw(positive negative);
    } else {
        my %recognized_ori = (
            plus => { plus  => 'positive', minus => 'unsupported', },
            minus => { plus  => 'negative', minus => 'unsupported', },
        );
        $ori = $recognized_ori{$one}{$two} || croak "Orientation $one - $two is not recognized";
        croak "Orientation $one - $two is $ori" unless grep {$ori eq $_} qw(positive negative);
    }
    
    return $ori;
}

# gets a list of asn "segs" to process, since some .asn files have nested or "disc-segs"
sub get_seg_list {
    my ($store, $type, $segs) = @_;

    if (ref($segs) eq 'HASH') {
        $segs = [$segs];
        croak "Expected type to be partial when there was only one seg, but type was $type" unless $type eq 'partial';
    } elsif (ref($segs) eq 'ARRAY') {
        $segs = [map {$_->{segs}} @{$segs}];
        croak "Expected type to be disc when there were many segs, but type was $type " unless $type eq 'disc';
    } else {
        croak "Not sure what kinds of segments are being processing: " . ref($segs);
    }

    my @ids = @{$segs->[0]{ids}};

    for my $subsegs (@{$segs}) {
        croak "No support for discontinuous segs nested within discontinuous segs"
            if ref($subsegs) eq 'ARRAY'; 
        croak sprintf("Not all disc segs had the same ids: %s and %s", join(',', @ids), join(',', @{$subsegs->{ids}}))
            if (@ids != @{$subsegs->{ids}});
        croak "Only 2 seqs are currently supported"
            if ($subsegs->{dim} != 2);
        # TODO 3 seqs are too complicated for now
        # TODO maybe verify types and dims here too?
    }
    
    return $segs;
}

# verify, process, and write a single remapping
sub process_remap_chunk {
    my $args = shift;

    my $store          = $args->{store};
    my $ori           = $args->{ori},
    my @alt           = @{$args->{alt}};
    my @par           = @{$args->{par}};
    my @len           = @{$args->{len}};
    my @ops           = @{$args->{ops}};
    my $par_fasta_seq = $args->{par_fasta_seq};
    my $alt_fasta_seq = $args->{alt_fasta_seq};
    my $par_fa        = $args->{par_fa};
    my $alt_fa        = $args->{alt_fa};
    my $count         = $args->{count};
    my $chromosome    = $args->{chromosome};
    my $outremapfh    = $args->{outremapfh};
    my $outfastafh    = $args->{outfastafh};

    my $flank = $store->{flank};

    croak "Sequence was first in alt scaffold, but the first cigar op is not a match" if
        ($ops[0] ne 'M' and $count == 0);
    croak "Sequence was after the first in alt scaffold, but did not skip the deletion" if
        ($ops[0] eq 'D' and $count != 0);
    croak "Last cigar op was not a match" if
        ($ops[$#ops] ne 'M');

    # TODO this assumes that the parent/primary sequence is always in the "positive" orientation
    # might not be an issue...
    my $parent_seq_length = $par_fa->length($par_fasta_seq->{seqid});

    my $parent_start = $par[0]->{start} + 1 - $flank;
    my $preflank_len = $flank;
    if ($parent_start < 1) {
        $preflank_len -= (1 - $parent_start);
        $parent_start = 1;
    }

    my $parent_stop = $par[$#par]->{stop} + $flank;
    my $postflank_len = $flank;
    if ($parent_stop > $parent_seq_length) {
        $postflank_len -= ($parent_stop - $parent_seq_length);
        $parent_stop = $parent_seq_length;
    }


    my $full_seqid   = join '_', grep {$_} split '\|', $alt_fasta_seq->{seqid};
    my $seq_line     = sprintf '>%s_%s', $full_seqid, $count;
    my $remap_seq_line = sprintf '%s-%s|%s|%s', $seq_line, $chromosome, $parent_start, $parent_stop;


    my $cigar = '';

    for (0 .. $#ops) {
        my $curlen = $len[$_];
        if ($_ == 0 and $ops[0] eq 'M') {
            $curlen += $preflank_len;
        } 
        if ($_ == $#ops and $ops[$#ops] eq 'M') {
            $curlen += $postflank_len;
        }
        $cigar .= $curlen.$ops[$_];
    }

    say $outremapfh $remap_seq_line;
    say $outremapfh $cigar;

    # make sure this is a "contiguous" segment
    my $alt_start;
    my $alt_stop;

    if ($ori eq 'positive') {
        for my $i (0 .. $#alt-1) {
            unless ( ($alt[$i]->{stop}) == $alt[$i+1]->{start} ) {
                say Dumper(\@alt);
                croak "Noncontiguous sequence at alt index $i!";
            }
        }

        $alt_start = $alt[0]->{start} + 1;
        $alt_stop  = $alt[$#alt]->{stop};
    } elsif ($ori eq 'negative') {
        for my $i (0 .. $#alt-1) {
            $i = $#alt - ($i + 1);
            unless ( ($alt[$i]->{start}) == $alt[$i+1]->{stop} ) {
                say Dumper(\@alt);
                croak "Noncontiguous sequence at alt index $i!";
            }
        }

        $alt_start = $alt[$#alt]->{start} + 1;
        $alt_stop  = $alt[0]->{stop};
    } else {
        croak "Nonsensical orientation '$ori'";
    }

    my $alt_length = $alt_stop - $alt_start + 1;
    my $seq        = $alt_fa->seq($alt_fasta_seq->{seqid}, $alt_start => $alt_stop);
    my $seq_length = length $seq;
    
    $seq = rev_complement($store, $seq) if ($ori eq 'negative');

    croak "Expected a seq of length $alt_length, but fasta gave $seq_length instead." if $alt_length != $seq_length;
    
    my $preflank_seq = $preflank_len > 0
        ? $par_fa->seq($par_fasta_seq->{seqid}, $parent_start => $parent_start - 1 + $preflank_len)
        : '';
    my $postflank_seq = $postflank_len > 0
        ? $par_fa->seq($par_fasta_seq->{seqid}, $parent_stop - $postflank_len + 1 => $parent_stop)
        : '';

    my $flanked_seq        = $preflank_seq . $seq . $postflank_seq;
    my $flanked_seq_length = length $flanked_seq;

    my $curpos     = 0;
    my $linelength = 70;

    say $outfastafh $seq_line.$alt_fasta_seq->{desc};
    
    while ($curpos < $flanked_seq_length) {
        say $outfastafh substr $flanked_seq, $curpos, $linelength;
        $curpos += $linelength;
    }
}

# reverse complement reads
sub rev_complement {
    my ($store, $seq) = @_;

    my %complement = (
        'A' => 'T',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
    );

    return join('', reverse(map {exists($complement{$_}) ? $complement{$_} : $_ } split('', $seq)));
}

# grammar for parsing an asn file
sub parse_asn {
    my ($store, $asn) = @_;
    
    my $data = read_file($asn) or croak "Failed to read $asn";

    #$::RD_TRACE = 1;

    # TODO there may be a way to make this grammar more efficient

    my $grammar = q!
        startrule: "Seq-align" "::=" body
                { $item[3] }

        end: /^\Z/
        
        body: '{' item(s /,/) '}'
                { $item[2] }
        
        item: type 
                { $return = {
                        'key'   => 'type',
                        'value' => $item{type}
                    }
                }
            | dim
                { $return = {
                        'key'   => 'dim',
                        'value' => $item{dim}
                    }
                }
            | score
                { $return = {
                        'key'   => 'score',
                        'value' => $item{score}
                    }
                }
            | segs
                { $return = {
                        'key'   => 'segs',
                        'value' => $item{segs}
                    }
                }
            
        type: "type" word
                { $item[2] }
        
        
        dim: "dim" int
                { $item[2] }
                
                
        score: "score" '{' keyval(s /,/) '}'
                { $item[3] }
            
        keyval: '{' key ',' value '}'
                { $return = {
                        'key'   => $item{key},
                        'value' => $item{value}
                    }
                }
                
        key: "id" "str" str
                { $item[3] }
        
        value: "value" "int" int
                { $item[3] }
            | "value" "real" real
                { $item[3] }
        
        
        # order is important here, note how the 2nd production includes segtype instead of terminal "disc"
        segs: "segs" "disc" '{' body(s /,/) '}'
                { $item[4] }
            | "segs" segtype '{' dim ',' numseg ',' segids ',' starts ',' lens ',' strands '}'
                { $return = {
                        'type'    => $item{segtype},
                        'dim'     => $item{dim},
                        'numseg'  => $item{numseg},
                        'ids'     => $item{segids},
                        'starts'  => $item{starts},
                        'lens'    => $item{lens},
                        'strands' => $item{strands}
                    }
                }
            | "segs" segtype '{' dim ',' numseg ',' segids ',' starts ',' lens '}'
                { $return = {
                        'type'   => $item{segtype},
                        'dim'    => $item{dim},
                        'numseg' => $item{numseg},
                        'ids'    => $item{segids},
                        'starts' => $item{starts},
                        'lens'   => $item{lens}
                    }
                }
        
        segtype: word
        
        numseg: "numseg" int
                { $item[2] }
                
        segids: "ids" '{' segid(s /,/) '}'
                { $item[3] }
                
        segid: "gi" int
                { $item[2] }
        
        starts: "starts" '{' int(s /,/) '}'
                { $item[3] }
        
        lens: "lens" '{' int(s /,/) '}'
                { $item[3] }
        
        strands: "strands" '{' strand(s /,/) '}'
                { $item[3] }
        
        strand: "plus" | "minus"
        
        
        str: '"' word '"'
                { $item[2] }
        
        int: /-?[0-9]+/
        
        real: '{' int ',' int ',' int '}'
                { [ $item[2], $item[4], $item[6] ] }
        
        word: /[a-zA-Z_\-0-9]+/ | backslash
        
        backslash: '\\\\' # TODO wtf?
    !;

    my $parser = new Parse::RecDescent($grammar) or croak "Bad grammar!";
    my $seqalign = $parser->startrule($data);
    croak "Failed to parse $asn" unless $seqalign;
    
    return verify_and_transform($store, $seqalign);
}

# verifies and transforms the data structure returned from the asn parse grammar
sub verify_and_transform {
    my ($store, $seqalign_orig) = @_;
    
    my $seqalign = transform_to_hash($store, $seqalign_orig);

    if ($seqalign->{type} eq 'not-set') {
        warn colorize('yellow', "Warning, Seq-align type marked as not-set, assuming 'partial'.\n");
        $seqalign->{type} = 'partial';
    };

    croak "Unrecognized type '$seqalign->{type}'" unless scalar(grep {$seqalign->{type} eq $_} qw(partial disc));
    
    # process scores
    my $scores = transform_to_hash($store, $seqalign->{score}) if $seqalign->{score};
    for (keys %{$scores}) {
        # convert reals
        $scores->{$_} = $scores->{$_}[0] * ($scores->{$_}[1] ** $scores->{$_}[2]) if (ref($scores->{$_}) eq 'ARRAY');
    }
    $seqalign->{score} = $scores;

    # recurse with discontinuous segs
    if ($seqalign->{type} eq 'disc') {
        croak "Expected multiple discontinuous segs" unless ref $seqalign->{segs} eq 'ARRAY';

        my @disc_segs;
        for (@{$seqalign->{segs}}) {
            push @disc_segs, verify_and_transform($store, $_);
        }

        $seqalign->{segs} = \@disc_segs;
    } else {
        croak "Mismatched dims found" if ($seqalign->{dim} and $seqalign->{dim} != $seqalign->{segs}{dim});
        croak "Expected partial segs to point to hash ref" unless ref $seqalign->{segs} eq 'HASH';

        $seqalign->{segs} = transform_seqalign_segs($store, $seqalign->{segs});
    }
    return $seqalign;
}

# helper for processing segments from parsed asn file
sub transform_seqalign_segs {
    my ($store, $seqalign_segs) = @_;

    my $dim     = $seqalign_segs->{dim};
    my $numseg  = $seqalign_segs->{numseg};
    
    if (!defined $seqalign_segs->{strands}) {
        warn colorize('yellow', "Warning, no strands parsed. Assuming all strands are positive. This may be undesirable behavior.\n");
        my $strands;
        for (0 .. $dim*$numseg) {
            push @$strands, 'plus';
        }
        $seqalign_segs->{strands} = $strands;
    }
    
    my ($numids, $numstarts, $numlens, $numstrands) = map { scalar(@{$seqalign_segs->{$_}}) } qw(ids starts lens strands);
    
    croak "Counts don't match"
        if ($numids != $dim
            and $numlens != $numseg
            and ($numstarts / $dim) != $numseg
            and ($numstrands / $dim) != $numseg);
    
    croak "Unrecognized segtype '$seqalign_segs->{type}'" unless scalar(grep {$seqalign_segs->{type} eq $_} qw(denseg));
    
    # transform starts and strands into array of arrays
    my $starts_ref = $seqalign_segs->{starts};
    my $strands_ref = $seqalign_segs->{strands};
    
    my @starts; my @strands;
    push @starts, [ splice @{$starts_ref}, 0, $dim ] while @{$starts_ref};
    push @strands, [ splice @{$strands_ref}, 0, $dim ] while @{$strands_ref};
    
    $seqalign_segs->{starts} = \@starts;
    $seqalign_segs->{strands} = \@strands;
    
    return $seqalign_segs;
}

# helper that transforms a list of key-value pairs into a hash (checks that keys are unique)
sub transform_to_hash {
    my ($store, $hash_ref) = @_;

    # verify keys are unique
    my %hash;
    for (@{$hash_ref}) {
        $hash{$_->{key}}++;
    } 

    croak "Repeat keys" if scalar(grep {$_ != 1} values %hash);

    # transform array of kv pairs into a plain hash
    for (@{$hash_ref}) {
        $hash{$_->{key}} = $_->{value};
    }

    return \%hash;
}

# colorize strings
sub colorize {
    my ($color, $string) = @_;
    my %colorvals = (
        "Running"     => 33,
        "Unstartable" => 35,
        "Succeeded"   => 32,
        "Failed"      => 31,
        "Scheduled"   => 34,
        "scheduled"   => 35,
        "blue"        => 34,
        "green"       => 32,
        "purple"      => 35,
        "gray"        => 30,
        "grey"        => 30,
        "Abandoned"   => 30,
        "new"         => 34,
        "running"     => 33,
        "yellow"      => 33,
        "done"        => 32,
        "crashed"     => 31,
        "red"         => 31,
        "clear"       => 0,
    );
    
    if ($string) {
        return "[".$colorvals{$color}."m".$string."[0m";
    } else {
        return "[".$colorvals{$color}."m";
    }
}


