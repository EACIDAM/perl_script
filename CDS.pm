#--------CDS.pm file------------#

package Cds;
use strict;

#constructeur
sub new_1
{
	my($class, $id, $start, $end) = @_ ;
	my $this={};

	bless($this, $class);

	$this->{ID}=$id;
	$this->{START}=$start;
	$this->{END}=$end;
	return $this;
}

sub new_2
{
	my ($classe) = @_ ;
	my $this={};
	bless($this, $classe) ;
	return $this;
}

#getters
sub get_id
{
	my($this) = @_ ;
	return $this->{ID};
}

sub get_start
{
	my($this) = @_ ;
	return $this->{START};
}

sub get_end
{
	my($this) = @_ ;
	return $this->{END};
}

#setters
sub set_id
{
	my($this, $id) = @_ ;
	$this->{ID}=$id;
}

sub set_start
{
	my($this, $start) = @_ ;
	$this->{START}=$start;
}

sub set_end
{
	my($this, $end) = @_ ;
	$this->{END}=$end;
}

1;
