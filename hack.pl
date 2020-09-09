$string = "";
while (<>) {
    chomp;
    if (/\D\d+\D/) {
      $string =~ s/\s+\+\/\-\s+/+\/-/g;
      $string =~ s/\s+/\t/g;
      print "$string\n";
      $string = "";
    }
    $string .= $_;
    /\/$/ || ($string .= " ");
}
