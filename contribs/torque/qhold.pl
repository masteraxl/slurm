#! /usr/bin/perl -w
################################################################################
#
# qhold - places a hold on moab jobs in familar pbs format.
#
#                 Copyright (c) 2006 Cluster Resources, Inc.
#
################################################################################

use strict;
use FindBin;
use Getopt::Long 2.24 qw(:config no_ignore_case);
use lib "${FindBin::Bin}/../tools";
use Moab::Tools;    # Required before including config to set $homeDir
use autouse 'Pod::Usage' => qw(pod2usage);
use XML::LibXML;
BEGIN { require "config.moab.pl"; }
our ($logLevel, $mjobctl);

Main:
{
    logPrint("Command line arguments: @ARGV\n") if $logLevel;

    # Parse Command Line Arguments
    my ($help, $hold, $man);
    GetOptions(
        'help|?' => \$help,
        'h=s'      => \$hold,
        'man'    => \$man,
      )
      or pod2usage(2);

    # Display usage if necessary
    pod2usage(0) if $help;
    if ($man)
    {
        if ($< == 0)    # Cannot invoke perldoc as root
        {
            my $id = eval { getpwnam("nobody") };
            $id = eval { getpwnam("nouser") } unless defined $id;
            $id = -2                          unless defined $id;
            $<  = $id;
        }
        $> = $<;                         # Disengage setuid
        $ENV{PATH} = "/bin:/usr/bin";    # Untaint PATH
        delete @ENV{'IFS', 'CDPATH', 'ENV', 'BASH_ENV'};
        if ($0 =~ /^([-\/\w\.]+)$/) { $0 = $1; }    # Untaint $0
        else { die "Illegal characters were found in \$0 ($0)\n"; }
        pod2usage(-exitstatus => 0, -verbose => 2);
    }

    # Use sole remaining argument as jobIds
    my @jobIds;
    if (@ARGV)
    {
        @jobIds = @ARGV;
    }
    else
    {
        pod2usage(2);
    }

    # Build command
    my $cmd = "$mjobctl -h ";
    if ($hold)
    {
        if    ($hold eq "u") { $cmd .= "user "; }
        elsif ($hold eq "s") { $cmd .= "system "; }
        elsif ($hold eq "o") { $cmd .= "batch "; }
        elsif ($hold eq "n") { exit 0; }
    }
    else
    {
        $cmd .= "user ";
    }
    $cmd .= join ' ', @jobIds;
    logPrint("Invoking subcommand: $cmd\n") if $logLevel >= 2;

    my $output = `$cmd 2>&1`;
    my $rc     = $?;

    # Display and log error output
    if ($output =~ /invalid job specified \(([^\(\)]+)\)/)
    {
        logWarn("qhold: Unknown Job Id $1\n");
    }
    elsif ($rc)
    {
        logWarn($output);
    }
    else
    {
        logPrint("Subcommand output:>\n$output") if $logLevel >= 3;
    }

    exit $rc;
}

##############################################################################

__END__

=head1 NAME

B<qhold> - places a hold on jobs in a familiar pbs format

=head1 SYNOPSIS

B<qhold> [B<-h> B<u>|B<o>|B<s>|B<n>] I<job_id>...

=head1 DESCRIPTION

The B<qhold> command requests that a hold be placed on a job.  A job that is on hold is not eligible for execution.  There are three supported holds: USER, OTHER (also known as operator or batch), and SYSTEM.

If the B<-h> option is not specified, the USER hold will be applied to the specified jobs.

=head1 OPTIONS

=over 4

=item B<-h> I<hold_type>

Specifieds the types of holds to be placed on the job.

The I<hold_type> argument is a one of the characters "u", "o", "s" or "n".  The hold type associated with each letter is:

    B<u> - USER

    B<o> - OTHER

    B<s> - SYSTEM

    B<n> - None

=item B<-? | --help>

brief help message

=item B<--man>

full documentation

=back

=head1 EXIT STATUS

On success, B<qhold> will exit with a value of zero. On failure, B<qhold> will exit with a value greater than zero.

=cut

