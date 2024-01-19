

document.querySelector('title').textContent = 'Perl on Biowulf';




 .hl { background-color: #ffff99; }
 code {padding: 1px; background-color: #eeeeee; border: 1px solid #bbbbbb;}
 dt {font-weight: bold; margin-top: 5px;}
 dd {padding-left: 2px; border-left: 1px solid #bbbbbb;}
 .btt {border: 1px solid silver;
 background-color: white;
 padding: 5px;
 position: relative;
 margin: 5px 0px 10px 10px;
 float: right;
 top: -25px;
 left: 10px;
 }

Perl on Biowulf


|  |
| --- |
| 
Quick Links
[Installed modules](#modules)
[Personal installation](#personal)
[Documentation](#docs)
 |




Perl is available by default on the biowulf cluster. More recent versions are also
available via the [modules](https://hpc.nih.gov/apps/modules.html) system.
There are a [large number of perl modules available](#modules) for these
versions.



In addition, users can install modules into their own /home directory as well.




Documentation
[top](#top)
* [Perl home](https://www.perl.org/)
* [Perl documentation](https://perldoc.perl.org/perl)
* [cpanm](https://metacpan.org/dist/App-cpanminus/view/bin/cpanm)



 
$(document).ready(function() { $('#perl-modules-ROCKY8').DataTable({
 columnDefs: [
 // Left align the header content of column 1
 { className: "dt-left", targets: [ '\_all' ] }
 ]
}); });

Installed modules
[top](#top)
The following modules are available for perl:






| Module | default | 5.18 | 5.34 | 5.36 | 5.38 |
| --- | --- | --- | --- | --- | --- |
| Acme::Comment | 1.04 | 1.04 | 1.04 | 1.04 | 1.04 |
| Acme::Damn | n/a | 0.08 | n/a | 0.08 | n/a |
| Algorithm::Combinatorics | 0.27 | 0.27 | 0.27 | 0.27 | 0.27 |
| Algorithm::Diff | 1.201 | 1.201 | 1.201 | 1.201 | 1.201 |
| Algorithm::Diff::XS | 0.04 | 0.04 | 0.04 | 0.04 | 0.04 |
| Algorithm::Munkres | 0.08 | 0.08 | 0.08 | 0.08 | 0.08 |
| aliased | 0.34 | 0.34 | 0.34 | 0.34 | 0.34 |
| Alien::Build | 2.48 | 2.80 | 2.74 | 2.80 | 2.80 |
| Alien::Build::Plugin::Download::GitLab | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 |
| Alien::Libxml2 | 0.17 | 0.19 | 0.19 | 0.19 | 0.19 |
| Any::URI::Escape | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 |
| Apache::LogFormat::Compiler | 0.36 | 0.36 | 0.36 | 0.36 | 0.36 |
| Apache::Session | 1.94 | 1.94 | 1.94 | 1.94 | 1.94 |
| AppConfig | 1.71 | 1.71 | 1.71 | 1.71 | 1.71 |
| App::cpanminus | 1.7046 | 1.7047 | 1.7046 | 1.7046 | 1.7046 |
| Archive::Any::Lite | 0.11 | 0.11 | 0.11 | 0.11 | 0.11 |
| Archive::Extract | 0.88 | 0.88 | 0.88 | 0.88 | 0.88 |
| Archive::Peek | 0.37 | 0.37 | 0.37 | 0.37 | 0.37 |
| Archive::Tar | 2.40 | 3.02 | 2.40 | 3.02 | 3.02 |
| Archive::Zip | 1.68 | 1.68 | 1.68 | 1.68 | 1.68 |
| Array::Compare | 3.0.8 | 3.0.8 | 3.0.8 | 3.0.8 | 3.0.8 |
| Array::Diff | 0.09 | 0.09 | 0.09 | 0.09 | 0.09 |
| Array::Unique | 0.08 | 0.09 | 0.08 | 0.09 | 0.09 |
| Array::Utils | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 |
| Authen::SASL | 2.16 | 2.1700 | 2.16 | 2.16 | 2.16 |
| Authen::Simple | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 |
| Authen::Simple::Passwd | 0.6 | 0.6 | 0.6 | 0.6 | 0.6 |
| autodie | 2.36 | 2.36 | n/a | n/a | n/a |
| bareword::filehandles | 0.007 | 0.007 | 0.007 | 0.007 | 0.007 |
| B::COW | 0.004 | 0.007 | 0.007 | 0.007 | 0.007 |
| B::Debug | n/a | n/a | 1.26 | 1.26 | 1.26 |
| Beam::Emitter | n/a | 1.007 | n/a | 1.007 | 1.007 |
| Beam::Service | n/a | 0.001 | n/a | 0.001 | 0.001 |
| Beam::Wire | n/a | 1.025 | n/a | 1.025 | 1.025 |
| B::Hooks::AtRuntime | n/a | n/a | n/a | 8 | 8 |
| B::Hooks::EndOfScope | 0.26 | 0.26 | 0.26 | 0.26 | 0.26 |
| B::Hooks::OP::Annotation | n/a | n/a | n/a | 0.44 | 0.44 |
| B::Hooks::OP::Check | 0.22 | 0.22 | 0.22 | 0.22 | 0.22 |
| bignum | 0.65 | 0.66 | 0.66 | 0.66 | n/a |
| Bio | undef | undef | undef | undef | undef |
| Bio::Align::Graphics | undef | undef | n/a | undef | undef |
| Bio::AlignIO::stockholm | undef | undef | undef | undef | undef |
| Bio::ASN1::EntrezGene | 1.73 | 1.73 | 1.73 | 1.73 | 1.73 |
| Bio::Cluster | 1.7.3 | 1.7.3 | 1.7.3 | 1.7.3 | 1.7.3 |
| Bio::LITE::Taxonomy | 0.07 | 0.07 | 0.07 | 0.07 | 0.07 |
| Bio::LITE::Taxonomy::NCBI | 0.1 | 0.1 | 0.1 | 0.1 | 0.1 |
| Bio::NeXMLIO | n/a | n/a | undef | undef | n/a |
| BioPerl | 1.7.8 | 1.7.8 | 1.7.8 | 1.7.8 | 1.7.8 |
| Bio::Phylo | v2.0.1 | v2.0.1 | v2.0.1 | v2.0.1 | v2.0.1 |
| Bio::SamTools | n/a | undef | undef | undef | undef |
| Bio::SearchIO::blastxml | 1.70 | 1.70 | 1.70 | 1.70 | 1.70 |
| Bio::Variation | 1.7.5 | 1.7.5 | 1.7.5 | 1.7.5 | 1.7.5 |
| Bit::Manip | 1.04 | 1.04 | 1.04 | 1.04 | 1.04 |
| Bit::Vector | 7.4 | 7.4 | 7.4 | 7.4 | 7.4 |
| boolean | 0.46 | 0.46 | 0.46 | 0.46 | 0.46 |
| BSD::Resource | 1.2911 | n/a | 1.2911 | 1.2911 | n/a |
| B::Utils | 0.27 | 0.27 | 0.27 | 0.27 | 0.27 |
| Cache::Cache | 1.08 | 1.08 | 1.08 | 1.08 | 1.08 |
| Cache::LRU | 0.04 | 0.04 | 0.04 | 0.04 | 0.04 |
| Cache::Simple::TimedExpiry | 0.27 | 0.27 | 0.27 | 0.27 | 0.27 |
| Canary::Stability | 2013 | 2013 | 2013 | 2013 | 2013 |
| Capture::Tiny | 0.48 | 0.48 | 0.48 | 0.48 | 0.48 |
| Carp | 1.50 | 1.50 | n/a | n/a | n/a |
| Carp::Assert | 0.21 | 0.22 | 0.21 | 0.22 | 0.22 |
| Carp::Clan | 6.08 | 6.08 | 6.08 | 6.08 | 6.08 |
| CGI | 4.54 | 4.57 | 4.54 | 4.57 | 4.57 |
| CGI::Emulate::PSGI | 0.23 | 0.23 | 0.23 | 0.23 | 0.23 |
| CGI::PSGI | 0.15 | 0.15 | 0.15 | 0.15 | 0.15 |
| Chart::Gnuplot | 0.23 | 0.23 | 0.23 | 0.23 | 0.23 |
| Chart::Gnuplot::Pie | 0.04 | 0.04 | 0.04 | 0.04 | 0.04 |
| Class::Accessor | 0.51 | 0.51 | 0.51 | 0.51 | 0.51 |
| Class::Accessor::Chained | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 |
| Class::Accessor::Lite | 0.08 | 0.08 | 0.08 | 0.08 | 0.08 |
| Class::Container | 0.13 | 0.13 | 0.13 | 0.13 | 0.13 |
| Class::Data::Inheritable | 0.09 | 0.09 | 0.09 | 0.09 | 0.09 |
| Class::Inspector | 1.36 | 1.36 | 1.36 | 1.36 | 1.36 |
| Class::Load | 0.25 | 0.25 | 0.25 | 0.25 | 0.25 |
| Class::Load::XS | 0.10 | 0.10 | 0.10 | 0.10 | 0.10 |
| Class::Measure | 0.09 | 0.10 | 0.09 | 0.10 | 0.10 |
| Class::MethodMaker | 2.24 | 2.24 | 2.24 | 2.24 | 2.24 |
| Class::Method::Modifiers | 2.13 | 2.15 | 2.13 | 2.15 | 2.15 |
| Class::Mix | 0.006 | 0.006 | 0.006 | 0.006 | 0.006 |
| Class::ReturnValue | 0.55 | 0.55 | 0.55 | 0.55 | 0.55 |
| Class::Singleton | 1.6 | 1.6 | 1.6 | 1.6 | 1.6 |
| Class::Std | 0.013 | 0.013 | 0.013 | 0.013 | 0.013 |
| Class::Std::Fast | 0.0.8 | 0.0.8 | 0.0.8 | 0.0.8 | 0.0.8 |
| Class::Tiny | 1.008 | 1.008 | 1.008 | 1.008 | 1.008 |
| Class::XSAccessor | 1.19 | 1.19 | 1.19 | 1.19 | 1.19 |
| Clone | 0.45 | 0.46 | 0.46 | 0.46 | 0.46 |
| Clone::Choose | 0.010 | 0.010 | 0.010 | 0.010 | 0.010 |
| common::sense | 3.75 | 3.75 | 3.75 | 3.75 | 3.75 |
| Compress::Raw::Bzip2 | 2.206 | 2.206 | 2.201 | 2.206 | 2.206 |
| Compress::Raw::Lzma | 2.103 | 2.206 | 2.201 | 2.204 | 2.206 |
| Compress::Raw::Zlib | 2.206 | 2.206 | 2.202 | 2.206 | 2.206 |
| Config::Any | n/a | 0.33 | n/a | 0.33 | 0.33 |
| Config::General | 2.65 | 2.65 | 2.65 | 2.65 | 2.65 |
| Config::Identity | 0.0019 | 0.0019 | 0.0019 | 0.0019 | 0.0019 |
| Config::IniFiles | 3.000003 | 3.000003 | 3.000003 | 3.000003 | 3.000003 |
| Config::Simple | 4.58 | 4.58 | 4.58 | 4.58 | 4.58 |
| Const::Fast | 0.014 | 0.014 | 0.014 | 0.014 | 0.014 |
| Convert::ASN1 | 0.33 | 0.34 | 0.33 | 0.33 | 0.33 |
| Convert::Binary::C | 0.84 | 0.84 | 0.84 | 0.84 | 0.84 |
| Convert::BinHex | 1.125 | 1.125 | 1.125 | 1.125 | 1.125 |
| Convert::Color | 0.11 | 0.17 | 0.13 | 0.17 | 0.17 |
| Cookie::Baker | 0.11 | 0.11 | 0.11 | 0.11 | 0.11 |
| CPAN | 2.34 | 2.36 | 2.34 | 2.36 | 2.36 |
| CPAN::DistnameInfo | n/a | 0.12 | 0.12 | 0.12 | 0.12 |
| Cpanel::JSON::XS | 4.27 | 4.37 | 4.32 | 4.36 | 4.37 |
| CPAN::FindDependencies | 3.13 | 3.13 | 3.13 | 3.13 | 3.13 |
| CPAN::Meta | n/a | 2.150010 | n/a | n/a | n/a |
| CPAN::Meta::Check | n/a | 0.018 | 0.014 | 0.017 | 0.018 |
| CPAN::Meta::Requirements | n/a | 2.143 | n/a | n/a | n/a |
| CPAN::Meta::YAML | n/a | 0.018 | n/a | n/a | n/a |
| CPANPLUS | 0.9914 | 0.9914 | 0.9914 | 0.9914 | 0.9914 |
| CPAN::SQLite | 0.220 | 0.220 | 0.220 | 0.220 | 0.220 |
| Crypt::Eksblowfish | 0.009 | 0.009 | 0.009 | 0.009 | 0.009 |
| Crypt::OpenSSL::Bignum | 0.09 | 0.09 | 0.09 | 0.09 | 0.09 |
| Crypt::OpenSSL::Guess | 0.15 | 0.15 | 0.15 | 0.15 | 0.15 |
| Crypt::OpenSSL::Random | 0.15 | 0.15 | 0.15 | 0.15 | 0.15 |
| Crypt::OpenSSL::RSA | 0.32 | 0.33 | 0.33 | 0.33 | 0.33 |
| Crypt::OpenSSL::X509 | 1.913 | 1.915 | 1.914 | 1.914 | 1.915 |
| Crypt::PasswdMD5 | 1.41 | 1.42 | 1.42 | 1.42 | 1.42 |
| Crypt::RC4 | 2.02 | 2.02 | 2.02 | 2.02 | 2.02 |
| CryptX | 0.076 | 0.078 | 0.077 | 0.078 | 0.078 |
| Crypt::X509 | 0.54 | 0.55 | 0.54 | 0.55 | 0.55 |
| CSS::Squish | 0.10 | 0.10 | 0.10 | 0.10 | 0.10 |
| Cwd | 3.75 | 3.75 | n/a | n/a | n/a |
| Cwd::Guard | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 |
| Data | undef | undef | undef | undef | undef |
| Data::Binary | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 |
| Data::Compare | 1.27 | 1.29 | 1.27 | 1.29 | 1.29 |
| Data::DPath | n/a | 0.59 | n/a | 0.58 | 0.59 |
| Data::Dump | 1.25 | 1.25 | 1.25 | 1.25 | 1.25 |
| Data::Dumper | n/a | 2.183 | n/a | n/a | n/a |
| Data::Dump::Streamer | 2.40 | 2.42 | 2.40 | 2.42 | 2.42 |
| Data::GUID | 0.050 | 0.051 | 0.050 | 0.051 | 0.051 |
| Data::ICal | 0.24 | 0.24 | 0.24 | 0.24 | 0.24 |
| Data::OptList | 0.113 | 0.114 | 0.112 | 0.114 | 0.114 |
| Data::Perl | 0.002011 | 0.002011 | 0.002011 | 0.002011 | 0.002011 |
| Data::Printer | n/a | n/a | n/a | 1.001000 | 1.001001 |
| Data::Record | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 |
| Data::Section | n/a | 0.200008 | 0.200007 | 0.200008 | 0.200008 |
| Data::Section::Simple | 0.07 | 0.07 | 0.07 | 0.07 | 0.07 |
| Data::TemporaryBag | 0.09 | 0.09 | 0.09 | 0.09 | 0.09 |
| Data::UUID | 1.226 | 1.226 | 1.226 | 1.226 | 1.226 |
| Data::Visitor | 0.31 | 0.32 | 0.31 | 0.32 | 0.32 |
| Date::Calc | 6.4 | 6.4 | 6.4 | 6.4 | 6.4 |
| Date::Calc::XS | 6.4 | 6.4 | 6.4 | 6.4 | 6.4 |
| Date::Extract | 0.06 | 0.07 | 0.06 | 0.07 | 0.07 |
| Date::Manip | 6.86 | 6.92 | 6.90 | 6.92 | 6.92 |
| Date::Parse | 2.33 | 2.33 | 2.33 | 2.33 | 2.33 |
| Date::Simple | 3.03 | 3.03 | 3.03 | 3.03 | 3.03 |
| DateTime | 1.58 | 1.59 | 1.59 | 1.59 | 1.59 |
| DateTime::Format::DateParse | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 |
| DateTime::Format::Mail | 0.403 | 0.403 | 0.403 | 0.403 | 0.403 |
| DateTime::Format::Natural | 1.13 | 1.17 | 1.13 | 1.16 | 1.17 |
| DateTime::Format::Strptime | 1.79 | 1.79 | 1.79 | 1.79 | 1.79 |
| DateTime::Format::W3CDTF | 0.08 | 0.08 | 0.08 | 0.08 | 0.08 |
| DateTime::HiRes | n/a | 0.04 | n/a | 0.04 | 0.04 |
| DateTime::Locale | 1.35 | 1.39 | 1.37 | 1.39 | 1.39 |
| DateTime::TimeZone | 2.52 | 2.60 | 2.53 | 2.53 | 2.53 |
| DBD::mysql | 4.050 | 4.050 | 4.050 | 4.050 | 4.050 |
| DBD::Pg | 3.16.0 | 3.17.0 | 3.16.0 | 3.16.3 | 3.16.3 |
| DBD::SQLite | 1.70 | 1.72 | 1.72 | 1.72 | 1.72 |
| DBI | 1.643 | 1.643 | 1.643 | 1.643 | 1.643 |
| DBIx::Connector | 0.57 | 0.59 | 0.58 | 0.58 | 0.59 |
| DBIx::DBSchema | 0.45 | 0.47 | 0.47 | 0.47 | 0.47 |
| DBIx::SearchBuilder | 1.71 | 1.78 | 1.74 | 1.76 | 1.78 |
| Devel::CallChecker | 0.008 | 0.009 | 0.008 | 0.009 | 0.009 |
| Devel::Caller | 2.06 | 2.07 | 2.06 | 2.07 | 2.07 |
| Devel::CheckCompiler | 0.07 | 0.07 | 0.07 | 0.07 | 0.07 |
| Devel::CheckLib | 1.14 | 1.16 | 1.16 | 1.16 | 1.16 |
| Devel::CheckOS | 1.93 | 1.96 | 1.95 | 1.96 | 1.96 |
| Devel::Cycle | 1.12 | 1.12 | 1.12 | 1.12 | 1.12 |
| Devel::GlobalDestruction | 0.14 | 0.14 | 0.14 | 0.14 | 0.14 |
| Devel::Hide | 0.0015 | 0.0015 | 0.0015 | 0.0015 | 0.0015 |
| Devel::LexAlias | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 |
| Devel::OverloadInfo | 0.007 | 0.007 | 0.007 | 0.007 | 0.007 |
| Devel::PPPort | 3.68 | 3.68 | n/a | n/a | n/a |
| Devel::ptkdb | 1.1091 | n/a | 1.1091 | 1.1091 | n/a |
| Devel::Required | 0.16 | 0.16 | 0.16 | 0.16 | 0.16 |
| Devel::Size | n/a | 0.83 | 0.83 | 0.83 | 0.83 |
| Devel::StackTrace | 2.04 | 2.04 | 2.04 | 2.04 | 2.04 |
| Devel::StackTrace::AsHTML | 0.15 | 0.15 | 0.15 | 0.15 | 0.15 |
| Devel::Symdump | 2.18 | 2.18 | 2.18 | 2.18 | 2.18 |
| Devel::Trace | 0.12 | 0.12 | 0.12 | 0.12 | 0.12 |
| Digest::BubbleBabble | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 |
| Digest::GOST | 0.06 | 0.06 | 0.06 | 0.06 | 0.06 |
| Digest::HMAC | 1.04 | 1.04 | 1.04 | 1.04 | 1.04 |
| Digest::MD5 | n/a | 2.58 | n/a | n/a | n/a |
| Digest::Perl::MD5 | 1.9 | 1.9 | 1.9 | 1.9 | 1.9 |
| Digest::SHA | 6.04 | n/a | n/a | n/a | n/a |
| Digest::SHA1 | 2.13 | 2.13 | 2.13 | 2.13 | 2.13 |
| DIME::Tools | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 |
| Dist::CheckConflicts | 0.11 | 0.11 | 0.11 | 0.11 | 0.11 |
| DynaLoader::Functions | 0.003 | 0.004 | 0.003 | 0.004 | 0.004 |
| Email::Abstract | 3.009 | n/a | n/a | n/a | n/a |
| Email::Address | 1.912 | 1.913 | 1.912 | 1.913 | 1.913 |
| Email::Address::List | 0.06 | 0.06 | 0.06 | 0.06 | 0.06 |
| Email::Address::XS | 1.04 | 1.05 | 1.05 | 1.05 | 1.05 |
| Email::Date::Format | 1.005 | 1.008 | 1.006 | 1.008 | 1.008 |
| Email::MessageID | 1.406 | 1.408 | 1.407 | 1.408 | 1.408 |
| Email::MIME | 1.952 | 1.953 | 1.952 | 1.953 | 1.953 |
| Email::MIME::ContentType | 1.026 | 1.028 | 1.027 | 1.028 | 1.028 |
| Email::MIME::Encodings | 1.315 | 1.317 | 1.316 | 1.317 | 1.317 |
| Email::Sender | 2.600 | n/a | n/a | n/a | n/a |
| Email::Simple | 2.216 | 2.218 | 2.216 | 2.218 | 2.218 |
| Encode | 3.17 | 3.19 | n/a | n/a | n/a |
| Encode::Locale | 1.05 | 1.05 | 1.05 | 1.05 | 1.05 |
| Env::Path | 0.19 | 0.19 | 0.19 | 0.19 | 0.19 |
| Error | n/a | 0.17029 | 0.17029 | 0.17029 | 0.17029 |
| Eval::Closure | 0.14 | 0.14 | 0.14 | 0.14 | 0.14 |
| Excel::Writer::XLSX | 1.09 | 1.11 | 1.09 | 1.11 | 1.11 |
| Exception::Class | 1.45 | 1.45 | 1.45 | 1.45 | 1.45 |
| experimental | n/a | 0.031 | n/a | n/a | n/a |
| Exporter::Lite | 0.09 | 0.09 | 0.09 | 0.09 | 0.09 |
| Exporter::Tiny | 1.002002 | 1.006002 | 1.006000 | 1.006002 | 1.006002 |
| ExtUtils::CBuilder | 0.280236 | 0.280236 | n/a | n/a | n/a |
| ExtUtils::CChecker | n/a | 0.11 | n/a | 0.11 | 0.11 |
| ExtUtils::Config | 0.008 | 0.008 | 0.008 | 0.008 | 0.008 |
| ExtUtils::Constant | 0.25 | 0.25 | 0.25 | 0.25 | 0.25 |
| ExtUtils::CppGuess | 0.26 | 0.26 | 0.26 | 0.26 | 0.26 |
| ExtUtils::Depends | 0.8001 | 0.8001 | 0.8001 | 0.8001 | 0.8001 |
| ExtUtils::Helpers | 0.026 | 0.026 | 0.026 | 0.026 | 0.026 |
| ExtUtils::InstallPaths | 0.012 | 0.012 | 0.012 | 0.012 | 0.012 |
| ExtUtils::MakeMaker | 7.64 | 7.70 | 7.64 | 7.70 | 7.70 |
| ExtUtils::MakeMaker::CPANfile | 0.09 | 0.09 | 0.09 | 0.09 | 0.09 |
| ExtUtils::Manifest | 1.75 | n/a | n/a | n/a | n/a |
| ExtUtils::ParseXS | 3.44 | 3.44 | 3.44 | n/a | n/a |
| ExtUtils::PkgConfig | 1.16 | 1.16 | 1.16 | 1.16 | 1.16 |
| FCGI | 0.82 | 0.82 | 0.82 | 0.82 | 0.82 |
| FCGI::Client | 0.09 | 0.09 | 0.09 | 0.09 | 0.09 |
| FCGI::ProcManager | 0.28 | 0.28 | 0.28 | 0.28 | 0.28 |
| FFI::CheckLib | 0.28 | 0.31 | 0.31 | 0.31 | 0.31 |
| File::chdir | 0.1010 | 0.1011 | 0.1011 | 0.1011 | 0.1011 |
| File::Copy::Link | 0.06 | 0.06 | 0.06 | 0.06 | 0.06 |
| File::Copy::Recursive | 0.45 | 0.45 | 0.45 | 0.45 | 0.45 |
| File::Copy::Recursive::Reduced | 0.006 | 0.006 | 0.006 | 0.006 | 0.006 |
| File::Find::Object | 0.3.6 | 0.3.8 | 0.3.6 | 0.3.7 | 0.3.8 |
| File::Find::Rule | 0.34 | 0.34 | 0.34 | 0.34 | 0.34 |
| File::Find::Rule::DirectoryEmpty | 1.11 | 1.11 | 1.11 | 1.11 | 1.11 |
| FileHandle::Unget | 0.1634 | 0.1634 | 0.1634 | 0.1634 | 0.1634 |
| File::HomeDir | n/a | 1.006 | 1.006 | 1.006 | 1.006 |
| File::Listing | 6.15 | 6.16 | 6.15 | 6.15 | 6.16 |
| File::Path | 2.18 | 2.18 | n/a | n/a | n/a |
| File::pushd | n/a | 1.016 | 1.016 | 1.016 | 1.016 |
| File::Remove | 1.60 | 1.61 | 1.61 | 1.61 | 1.61 |
| File::ShareDir | 1.118 | 1.118 | 1.118 | 1.118 | 1.118 |
| File::ShareDir::Install | 0.13 | 0.14 | 0.14 | 0.14 | 0.14 |
| File::Slurp | 9999.32 | 9999.32 | 9999.32 | 9999.32 | 9999.32 |
| File::Slurper | 0.013 | 0.014 | 0.013 | 0.014 | 0.014 |
| File::Slurp::Tiny | 0.004 | 0.004 | 0.004 | 0.004 | 0.004 |
| File::Sort | 1.01 | 1.01 | 1.01 | 1.01 | 1.01 |
| File::Sync | 0.11 | 0.11 | 0.11 | 0.11 | 0.11 |
| Filesys::Notify::Simple | 0.14 | 0.14 | 0.14 | 0.14 | 0.14 |
| File::Tail | 1.3 | 1.3 | 1.3 | 1.3 | 1.3 |
| File::Tee | 0.07 | 0.07 | 0.07 | 0.07 | 0.07 |
| File::Temp | 0.2311 | n/a | n/a | n/a | n/a |
| File::Touch | 0.12 | 0.12 | 0.12 | 0.12 | 0.12 |
| File::TreeCreate | 0.0.1 | 0.0.1 | 0.0.1 | 0.0.1 | 0.0.1 |
| File::Type | 0.22 | 0.22 | 0.22 | 0.22 | 0.22 |
| File::Which | n/a | 1.27 | 1.27 | 1.27 | 1.27 |
| Font::AFM | 1.20 | 1.20 | 1.20 | 1.20 | 1.20 |
| Font::TTF | 1.06 | 1.06 | 1.06 | 1.06 | 1.06 |
| forks | n/a | 0.36 | n/a | 0.36 | n/a |
| Function::Parameters | n/a | n/a | n/a | 2.002003 | 2.002004 |
| GD | 2.76 | 2.78 | 2.76 | 2.77 | 2.78 |
| GD::Graph | 1.54 | 1.56 | 1.54 | 1.56 | 1.56 |
| GD::Graph::histogram | 1.1 | 1.1 | 1.1 | 1.1 | 1.1 |
| GD::Text | 0.86 | 0.86 | 0.86 | 0.86 | 0.86 |
| Geo::Distance | 0.25 | 0.25 | 0.25 | 0.25 | 0.25 |
| Getopt::Long | n/a | 2.54 | n/a | n/a | n/a |
| Getopt::Long::Descriptive | 0.110 | 0.111 | 0.110 | 0.111 | 0.111 |
| Getopt::Simple | 1.52 | 1.52 | 1.52 | 1.52 | 1.52 |
| GIS::Distance | 0.19 | 0.20 | 0.19 | 0.20 | 0.20 |
| GIS::Distance::Fast | 0.15 | 0.16 | 0.15 | 0.16 | 0.16 |
| GnuPG::Interface | 1.02 | 1.02 | 1.02 | 1.02 | 1.02 |
| GO | undef | undef | undef | undef | undef |
| Graph | 0.9725 | 0.9704 | 0.9704 | 0.9704 | 0.9704 |
| Graphics::ColorUtils | 0.17 | 0.17 | 0.17 | 0.17 | 0.17 |
| GraphViz | 2.24 | 2.26 | 2.25 | 2.26 | 2.26 |
| GSSAPI | 0.28 | 0.28 | 0.28 | 0.28 | 0.28 |
| Hash::FieldHash | 0.15 | 0.15 | 0.15 | 0.15 | 0.15 |
| Hash::Merge | 0.302 | n/a | 0.302 | 0.302 | 0.302 |
| Hash::MultiValue | 0.16 | 0.16 | 0.16 | 0.16 | 0.16 |
| Heap | 0.80 | 0.80 | 0.80 | 0.80 | 0.80 |
| Hijk | 0.28 | 0.28 | 0.28 | 0.28 | 0.28 |
| Hook::LexWrap | 0.26 | 0.26 | 0.26 | 0.26 | 0.26 |
| HTML::Formatter | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 |
| HTML::FormatText::WithLinks | 0.15 | 0.15 | 0.15 | 0.15 | 0.15 |
| HTML::FormatText::WithLinks::AndTables | 0.07 | 0.07 | 0.07 | 0.07 | 0.07 |
| HTML::Mason | 1.59 | 1.60 | 1.59 | 1.60 | 1.60 |
| HTML::Mason::PSGIHandler | 0.53 | 0.53 | 0.53 | 0.53 | 0.53 |
| HTML::Parser | 3.78 | 3.81 | 3.80 | 3.81 | 3.81 |
| HTML::Quoted | 0.04 | 0.05 | 0.04 | 0.04 | 0.05 |
| HTML::RewriteAttributes | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 |
| HTML::Scrubber | 0.19 | 0.19 | 0.19 | 0.19 | 0.19 |
| HTML::Table | undef | undef | undef | undef | undef |
| HTML-TableExtract | undef | undef | undef | undef | undef |
| HTML::Tagset | 3.20 | 3.20 | 3.20 | 3.20 | 3.20 |
| HTML::Template | 2.97 | 2.97 | 2.97 | 2.97 | 2.97 |
| HTML::Tree | 5.07 | 5.07 | 5.07 | 5.07 | 5.07 |
| HTTP::CookieJar | 0.012 | 0.014 | 0.014 | 0.014 | 0.014 |
| HTTP::Cookies | 6.10 | 6.10 | 6.10 | 6.10 | 6.10 |
| HTTP::Daemon | 6.14 | 6.16 | 6.14 | 6.16 | 6.16 |
| HTTP::Date | 6.05 | 6.06 | 6.05 | 6.05 | 6.06 |
| HTTP::Entity::Parser | 0.25 | 0.25 | 0.25 | 0.25 | 0.25 |
| HTTP::Headers::Fast | 0.22 | 0.22 | 0.22 | 0.22 | 0.22 |
| HTTP::Message | 6.36 | 6.44 | 6.44 | 6.44 | 6.44 |
| HTTP::MultiPartParser | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 |
| HTTP::Negotiate | 6.01 | 6.01 | 6.01 | 6.01 | 6.01 |
| HTTP::Request::AsCGI | 1.2 | 1.2 | 1.2 | 1.2 | 1.2 |
| HTTP::Server::Simple | 0.52 | 0.52 | 0.52 | 0.52 | 0.52 |
| HTTP::Server::Simple::PSGI | 0.16 | 0.16 | 0.16 | 0.16 | 0.16 |
| HTTP::Tiny | 0.080 | 0.088 | n/a | n/a | n/a |
| Image::Size | 3.300 | 3.300 | 3.300 | 3.300 | 3.300 |
| Importer | 0.026 | 0.026 | 0.026 | 0.026 | 0.026 |
| Import::Into | n/a | n/a | n/a | 1.002005 | 1.002005 |
| indirect | 0.39 | 0.39 | 0.39 | 0.39 | 0.39 |
| Inline | 0.86 | 0.86 | 0.86 | 0.86 | 0.86 |
| Inline::C | 0.82 | 0.82 | 0.82 | 0.82 | 0.82 |
| install | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 |
| Interpolation | 0.74 | 0.74 | 0.74 | 0.74 | 0.74 |
| IO::Capture | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 |
| IO::CaptureOutput | 1.1105 | 1.1105 | 1.1105 | 1.1105 | 1.1105 |
| IO::Compress | undef | undef | undef | undef | undef |
| IO::Compress::Lzma | 2.103 | 2.206 | 2.201 | 2.204 | 2.206 |
| IO::HTML | 1.004 | 1.004 | 1.004 | 1.004 | 1.004 |
| IO::SessionData | 1.03 | 1.03 | 1.03 | 1.03 | 1.03 |
| IO::Socket::INET6 | 2.73 | 2.73 | 2.73 | 2.73 | 2.73 |
| IO::Socket::IP | n/a | 0.42 | n/a | n/a | n/a |
| IO::Socket::SSL | n/a | 2.083 | 2.078 | 2.083 | 2.083 |
| IO::Socket::Timeout | 0.32 | 0.32 | 0.32 | 0.32 | 0.32 |
| IO::String | 1.08 | 1.08 | 1.08 | 1.08 | 1.08 |
| IO::Stringy | 2.113 | 2.113 | 2.113 | 2.113 | 2.113 |
| IO::Tty | 1.16 | 1.17 | 1.17 | 1.17 | 1.17 |
| IO::Zlib | 1.14 | n/a | n/a | 1.14 | n/a |
| IPC::Run | 20200505.0 | 20220807.0 | 20220807.0 | 20220807.0 | 20220807.0 |
| IPC::Run3 | 0.048 | 0.048 | 0.048 | 0.048 | 0.048 |
| IPC::Signal | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| IPC::System::Simple | n/a | 1.30 | 1.30 | 1.30 | 1.30 |
| IRI | 0.011 | 0.011 | 0.011 | 0.011 | 0.011 |
| Iterator | n/a | 0.03 | n/a | 0.03 | 0.03 |
| Iterator::Util | n/a | 0.02 | n/a | 0.02 | 0.02 |
| JSON | 4.05 | 4.10 | 4.10 | 4.10 | 4.10 |
| JSON::Any | 1.39 | 1.39 | 1.39 | 1.39 | 1.39 |
| JSON::MaybeXS | 1.004003 | 1.004005 | 1.004004 | 1.004005 | 1.004005 |
| JSON::PP | 4.16 | 4.16 | n/a | n/a | n/a |
| JSON::XS | 4.03 | 4.03 | 4.03 | 4.03 | 4.03 |
| Lexical::SealRequireHints | 0.011 | 0.012 | 0.011 | 0.012 | 0.012 |
| libintl-perl | undef | undef | undef | undef | undef |
| libwww::perl | undef | n/a | n/a | n/a | n/a |
| libxml-enno | n/a | n/a | undef | n/a | n/a |
| libxml-perl | undef | undef | undef | undef | undef |
| Lingua::EN::Inflect | 1.905 | 1.905 | 1.905 | 1.905 | 1.905 |
| Lingua::EN::Sentence | 0.31 | 0.34 | 0.33 | 0.33 | 0.34 |
| Lingua::Sentence | 1.100 | 1.100 | 1.100 | 1.100 | 1.100 |
| Linux::ACL | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 |
| Linux::Inotify2 | 2.3 | 2.3 | 2.3 | 2.3 | 2.3 |
| List::Compare | 0.55 | 0.55 | 0.55 | 0.55 | 0.55 |
| List::MoreUtils | 0.430 | 0.430 | 0.430 | 0.430 | 0.430 |
| List::MoreUtils::XS | 0.430 | 0.430 | 0.430 | 0.430 | 0.430 |
| List::Util | 1.62 | 1.63 | 1.63 | n/a | n/a |
| List::UtilsBy | 0.11 | 0.12 | 0.12 | 0.12 | 0.12 |
| load | 0.25 | 0.25 | 0.25 | 0.25 | 0.25 |
| Locale::Codes | 3.70 | n/a | 3.72 | 3.74 | 3.74 |
| Locale::Maketext::Fuzzy | 0.11 | 0.11 | 0.11 | 0.11 | 0.11 |
| Locale::Maketext::Lexicon | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| Log::Any | 1.710 | 1.717 | 1.713 | 1.715 | 1.716 |
| Log::Any::Adapter::Callback | 0.101 | 0.101 | 0.101 | 0.101 | 0.101 |
| Log::Dispatch | 2.70 | 2.71 | 2.70 | 2.71 | 2.71 |
| Log::Dispatch::Array | 1.003 | 1.005 | 1.003 | 1.005 | 1.005 |
| Logger::Simple | 2.0 | 2.0 | 2.0 | 2.0 | 2.0 |
| Log::Log4perl | 1.54 | 1.57 | 1.57 | 1.57 | 1.57 |
| Log::Message | 0.08 | 0.08 | 0.08 | 0.08 | 0.08 |
| Log::Message::Simple | 0.10 | 0.10 | 0.10 | 0.10 | 0.10 |
| Log::Report | 1.33 | 1.34 | 1.34 | 1.34 | 1.34 |
| Log::Report::Optional | 1.07 | 1.07 | 1.07 | 1.07 | 1.07 |
| LWP | 6.71 | 6.72 | 6.71 | 6.71 | 6.72 |
| LWP::MediaTypes | 6.04 | 6.04 | 6.04 | 6.04 | 6.04 |
| LWP::Protocol::http10 | 6.03 | 6.03 | 6.03 | 6.03 | 6.03 |
| LWP::Protocol::https | 6.10 | 6.11 | 6.10 | 6.10 | 6.11 |
| Mail::Mbox::MessageParser | 1.5111 | 1.5111 | 1.5111 | 1.5111 | 1.5111 |
| MailTools | 2.21 | 2.21 | 2.21 | 2.21 | 2.21 |
| Math::Bezier | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 |
| Math::Big | 1.16 | 1.16 | 1.16 | 1.16 | 1.16 |
| Math::BigInt | 1.999830 | 1.999839 | 1.999837 | 1.999842 | 1.999842 |
| Math::BigRat | 0.2622 | 0.2624 | 0.2624 | 0.2624 | n/a |
| Math::CDF | 0.1 | 0.1 | 0.1 | 0.1 | 0.1 |
| Math::Derivative | 1.01 | 1.01 | 1.01 | 1.01 | 1.01 |
| Math::Random | 0.72 | 0.72 | 0.72 | 0.72 | 0.72 |
| Math::Random::MT::Auto | 6.23 | 6.23 | 6.23 | 6.23 | 6.23 |
| Math::Round | 0.07 | 0.07 | 0.07 | 0.07 | 0.07 |
| Math::Spline | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 |
| Math::Utils | 1.14 | 1.14 | 1.14 | 1.14 | 1.14 |
| Math::VecStat | 0.08 | 0.08 | 0.08 | 0.08 | 0.08 |
| MCE | 1.878 | 1.888 | 1.882 | 1.887 | 1.888 |
| Memoize | n/a | n/a | n/a | n/a | 1.03 |
| Memory::Usage | 0.201 | 0.201 | 0.201 | 0.201 | 0.201 |
| MIME::Charset | 1.012.2 | 1.013.1 | 1.013.1 | 1.013.1 | 1.013.1 |
| MIME::Lite | 3.033 | 3.033 | 3.033 | 3.033 | 3.033 |
| MIME::Lite::TT::HTML | 0.04 | 0.04 | 0.04 | 0.04 | 0.04 |
| MIME::Tools | 5.509 | 5.510 | 5.510 | 5.510 | 5.510 |
| MIME::Types | 2.22 | 2.24 | 2.22 | 2.24 | 2.24 |
| Mixin::Linewise | 0.110 | 0.111 | 0.110 | 0.111 | 0.111 |
| Mock::Config | 0.03 | 0.03 | 0.03 | 0.03 | 0.03 |
| Module::Build | 0.4231 | 0.4234 | 0.4232 | 0.4234 | 0.4234 |
| Module::Build::Tiny | 0.039 | 0.046 | 0.039 | 0.046 | 0.046 |
| Module::Build::XSUtil | 0.19 | 0.19 | 0.19 | 0.19 | 0.19 |
| Module::CoreList | 5.20221120 | 5.20230820 | 5.20180120 | 5.20230423 | 5.20230720 |
| Module::CPANfile | 1.1004 | 1.1004 | 1.1004 | 1.1004 | 1.1004 |
| Module::CPANTS::Analyse | 1.01 | 1.01 | 1.01 | 1.01 | 1.01 |
| Module::ExtractUse | 0.344 | 0.345 | 0.344 | 0.345 | 0.345 |
| Module::Find | 0.15 | 0.16 | 0.16 | 0.16 | 0.16 |
| Module::Implementation | 0.09 | 0.09 | 0.09 | 0.09 | 0.09 |
| Module::Install | 1.19 | 1.21 | 1.19 | 1.21 | 1.21 |
| Module::Install::CPANfile | 0.12 | 0.12 | 0.12 | 0.12 | 0.12 |
| Module::Load | 0.36 | 0.36 | n/a | n/a | n/a |
| Module::Load::Conditional | n/a | 0.74 | n/a | n/a | n/a |
| Module::Metadata | n/a | 1.000038 | n/a | n/a | n/a |
| Module::Pluggable | 5.2 | 5.2 | 5.2 | 5.2 | 5.2 |
| Module::Refresh | 0.17 | 0.18 | 0.18 | 0.18 | 0.18 |
| Module::Runtime | 0.016 | 0.016 | 0.016 | 0.016 | 0.016 |
| Module::Runtime::Conflicts | 0.003 | 0.003 | 0.003 | 0.003 | 0.003 |
| Module::ScanDeps | 1.31 | 1.33 | 1.31 | 1.31 | 1.32 |
| Module::Signature | 0.88 | 0.88 | 0.88 | 0.88 | 0.88 |
| Module::Util | 1.09 | 1.09 | 1.09 | 1.09 | 1.09 |
| Module-Versions-Report | undef | undef | undef | undef | undef |
| Mojo::IOLoop::Delay | 8.76 | n/a | undef | n/a | n/a |
| Mojolicious | 9.24 | 9.33 | 7.94 | 9.32 | 9.33 |
| Moo | 2.004004 | 2.004004 | 2.004004 | 2.004004 | 2.004004 |
| Moose | 2.2201 | 2.2206 | 2.2201 | 2.2203 | 2.2206 |
| MooseX::ArrayRef | 0.005 | 0.005 | 0.005 | 0.005 | 0.005 |
| MooseX::Extended | n/a | n/a | n/a | 0.35 | 0.35 |
| MooseX::Role::Parameterized | 1.11 | 1.11 | 1.11 | 1.11 | 1.11 |
| MooseX::Role::WarnOnConflict | n/a | n/a | n/a | 0.01 | 0.01 |
| MooseX::Role::WithOverloading | 0.17 | 0.17 | 0.17 | 0.17 | 0.17 |
| MooseX::StrictConstructor | n/a | n/a | n/a | 0.21 | 0.21 |
| MooseX::Types | 0.50 | 0.50 | 0.50 | 0.50 | 0.50 |
| MooseX::Types::Path::Class | 0.09 | 0.09 | 0.09 | 0.09 | 0.09 |
| MooX::HandlesVia | 0.001009 | 0.001009 | 0.001009 | 0.001009 | 0.001009 |
| MooX::late | 0.100 | 0.100 | 0.100 | 0.100 | 0.100 |
| MooX::Locale::Passthrough | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 |
| MooX::Options | 4.103 | 4.103 | 4.103 | 4.103 | 4.103 |
| MooX::ProtectedAttributes | n/a | 0.03 | n/a | 0.03 | 0.03 |
| MooX::Types::MooseLike | 0.29 | 0.29 | 0.29 | 0.29 | 0.29 |
| MooX::TypeTiny | 0.002003 | 0.002003 | 0.002003 | 0.002003 | 0.002003 |
| Mouse | v2.5.10 | v2.5.10 | v2.5.10 | v2.5.10 | v2.5.10 |
| Mozilla::CA | 20211001 | 20230821 | 20221114 | 20221114 | 20221114 |
| Mozilla::PublicSuffix | v1.0.6 | v1.0.6 | v1.0.6 | v1.0.6 | v1.0.6 |
| MRO::Compat | 0.15 | 0.15 | 0.15 | 0.15 | 0.15 |
| multidimensional | 0.014 | 0.014 | 0.014 | 0.014 | 0.014 |
| MySQL::Config | 1.04 | 1.04 | 1.04 | 1.04 | 1.04 |
| namespace::autoclean | 0.29 | 0.29 | 0.29 | 0.29 | 0.29 |
| namespace::clean | 0.27 | 0.27 | 0.27 | 0.27 | 0.27 |
| Net | undef | undef | undef | undef | n/a |
| Net::CIDR | 0.21 | 0.21 | 0.21 | 0.21 | 0.21 |
| Net::Curl | 0.50 | 0.54 | 0.52 | 0.53 | 0.54 |
| Net::DNS | 1.33 | 1.20 | 1.20 | 1.20 | 1.20 |
| Net::HTTP | 6.22 | 6.23 | 6.22 | 6.22 | 6.23 |
| Net::IP | 1.26 | 1.26 | 1.26 | 1.26 | 1.26 |
| Net::LDAP | 0.68 | 0.68 | 0.68 | 0.68 | 0.68 |
| Net::OAuth | 0.28 | 0.28 | 0.28 | 0.28 | 0.28 |
| Net::Ping | 2.75 | n/a | n/a | n/a | n/a |
| Net::Server | 2.010 | 2.014 | 2.013 | 2.014 | 2.014 |
| Net::SMTP::TLS | 0.12 | 0.12 | 0.12 | 0.12 | 0.12 |
| Net::SSLeay | n/a | 1.92 | 1.92 | 1.92 | 1.92 |
| Net::Twitter | 4.01043 | 4.01043 | 4.01043 | 4.01043 | 4.01043 |
| Number::Compare | 0.03 | 0.03 | 0.03 | 0.03 | 0.03 |
| Number::Format | 1.75 | 1.76 | 1.75 | 1.76 | 1.76 |
| Number::Range | 0.12 | 0.12 | 0.12 | 0.12 | 0.12 |
| Object::Accessor | 0.48 | 0.48 | 0.48 | 0.48 | 0.48 |
| Object::InsideOut | 4.05 | 4.05 | 4.05 | 4.05 | 4.05 |
| Object::Pad | n/a | n/a | n/a | 0.79 | 0.79 |
| OLE::Storage\_Lite | 0.20 | 0.22 | 0.20 | 0.22 | 0.22 |
| Package::Constants | 0.06 | n/a | 0.06 | 0.06 | 0.06 |
| Package::DeprecationManager | 0.17 | 0.18 | 0.17 | 0.18 | 0.18 |
| Package::Stash | 0.40 | 0.40 | 0.40 | 0.40 | 0.40 |
| Package::Stash::XS | 0.29 | 0.30 | 0.30 | 0.30 | 0.30 |
| PadWalker | 2.5 | 2.5 | 2.5 | 2.5 | 2.5 |
| Parallel::ForkManager | 2.02 | 2.02 | 2.02 | 2.02 | 2.02 |
| Parallel::Prefork | 0.18 | 0.18 | 0.18 | 0.18 | 0.18 |
| Params::Classify | 0.015 | 0.015 | 0.015 | 0.015 | 0.015 |
| Params::Util | 1.102 | 1.102 | 1.102 | 1.102 | 1.102 |
| Params::Validate | 1.30 | 1.31 | 1.31 | 1.31 | 1.31 |
| Params::ValidationCompiler | 0.30 | 0.31 | 0.30 | 0.31 | 0.31 |
| PAR::Dist | 0.51 | 0.52 | 0.51 | 0.52 | 0.52 |
| Parse::CPAN::Packages | 2.40 | 2.40 | 2.40 | 2.40 | 2.40 |
| Parse::Distname | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 |
| Parse::LocalDistribution | 0.19 | 0.19 | 0.19 | 0.19 | 0.19 |
| Parse::PMFile | n/a | 0.44 | 0.43 | 0.44 | 0.44 |
| Parse::RecDescent | 1.967015 | 1.967015 | 1.967015 | 1.967015 | 1.967015 |
| Parse::Yapp | 1.21 | 1.21 | 1.21 | 1.21 | 1.21 |
| Path::Class | 0.37 | 0.37 | 0.37 | 0.37 | 0.37 |
| Path::Tiny | 0.122 | 0.144 | 0.144 | 0.144 | 0.144 |
| PAUSE::Permissions | 0.17 | 0.17 | 0.17 | 0.17 | 0.17 |
| PAUSE::Permissions::MetaCPAN | 0.001 | 0.100 | 0.001 | 0.100 | 0.100 |
| PDF::API2 | 2.043 | 2.044 | 2.044 | 2.044 | 2.044 |
| Pegex | 0.75 | 0.75 | 0.75 | 0.75 | 0.75 |
| Perl | 5.26.3 | 5.18.4 | 5.34.1 | 5.36.1 | 5.38.0 |
| PerlIO::eol | 0.17 | 0.19 | 0.17 | 0.18 | 0.19 |
| PerlIO::gzip | 0.20 | 0.20 | 0.20 | 0.20 | 0.20 |
| PerlIO::utf8\_strict | 0.009 | 0.010 | 0.010 | 0.010 | 0.010 |
| PerlIO::via::Timeout | 0.32 | 0.32 | 0.32 | 0.32 | 0.32 |
| Perl::PrereqScanner::NotQuiteLite | 0.9916 | 0.9917 | 0.9916 | 0.9917 | 0.9917 |
| Perl::Tidy | 20220217 | 20230701 | 20221112 | 20230309 | 20230701 |
| Perl::Unsafe::Signals | n/a | n/a | n/a | 0.03 | n/a |
| Plack | 1.0048 | 1.0050 | 1.0050 | 1.0050 | 1.0050 |
| Pod::Coverage | 0.23 | 0.23 | 0.23 | 0.23 | 0.23 |
| Pod::Coverage::TrustPod | 0.100005 | 0.100006 | 0.100005 | 0.100006 | 0.100006 |
| Pod::Escapes | n/a | 1.07 | n/a | n/a | n/a |
| Pod::Eventual | 0.094002 | 0.094003 | 0.094002 | 0.094003 | 0.094003 |
| Pod::LaTeX | 0.61 | n/a | 0.61 | 0.61 | 0.61 |
| Pod::Parser | n/a | 1.66 | 1.65 | 1.66 | 1.66 |
| Pod::Perldoc | n/a | 3.28 | n/a | n/a | n/a |
| Pod::Simple | 3.43 | 3.45 | 3.43 | 3.45 | 3.45 |
| Pod::Spell | 1.20 | 1.26 | 1.25 | 1.26 | 1.26 |
| Pod::Spell::CommonMistakes | 1.002 | 1.002 | 1.002 | 1.002 | 1.002 |
| Pod::Strip | 1.100 | 1.100 | 1.100 | 1.100 | 1.100 |
| POSIX::strftime::Compiler | 0.44 | 0.44 | 0.44 | 0.44 | 0.44 |
| PostScript::Metrics | 0.06 | 0.06 | 0.06 | 0.06 | 0.06 |
| PPI | 1.273 | 1.276 | 1.276 | 1.276 | 1.276 |
| Proc::Background | 1.30 | 1.32 | 1.30 | 1.32 | 1.32 |
| Proc::ProcessTable | 0.634 | 0.636 | 0.634 | 0.635 | 0.636 |
| Proc::Wait3 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 |
| Readonly | 2.05 | 2.05 | 2.05 | 2.05 | 2.05 |
| Redis | 1.999 | 2.000 | 1.999 | 2.000 | 2.000 |
| Ref::Util | 0.204 | 0.204 | 0.204 | 0.204 | 0.204 |
| Ref::Util::XS | 0.117 | 0.117 | 0.117 | 0.117 | 0.117 |
| Regexp::Common | 2017060201 | 2017060201 | 2017060201 | 2017060201 | 2017060201 |
| Regexp::Common::net::CIDR | 0.03 | 0.03 | 0.03 | 0.03 | 0.03 |
| Regexp::IPv6 | 0.03 | 0.03 | 0.03 | 0.03 | 0.03 |
| Regexp::Trie | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 |
| Regexp::Util | 0.005 | 0.005 | 0.005 | 0.005 | 0.005 |
| RL | n/a | 0.09 | n/a | 0.09 | 0.09 |
| Role::Basic | 0.13 | 0.13 | 0.13 | 0.13 | 0.13 |
| Role::Hooks | n/a | 0.008 | 0.008 | 0.008 | 0.008 |
| Role::Tiny | 2.002004 | 2.002004 | 2.002004 | 2.002004 | 2.002004 |
| RTF::Writer | 1.11 | 1.11 | 1.11 | 1.11 | 1.11 |
| Scope::Guard | 0.21 | 0.21 | 0.21 | 0.21 | 0.21 |
| Search::Elasticsearch | 6.81 | 6.81 | 6.81 | 6.81 | 6.81 |
| Search::Elasticsearch::Client::6\_0 | n/a | 6.81 | 6.81 | 6.81 | 6.81 |
| Search::Elasticsearch::Client::7\_0 | n/a | 7.713 | n/a | 7.713 | 7.713 |
| Sereal::Decoder | 4.023 | 5.004 | 5.001 | 5.004 | 5.004 |
| Sereal::Encoder | 4.023 | 5.004 | 5.001 | 5.004 | 5.004 |
| Server::Starter | 0.35 | 0.35 | 0.35 | 0.35 | 0.35 |
| Set::IntervalTree | 0.12 | 0.12 | 0.12 | 0.12 | 0.12 |
| Set::IntSpan | 1.19 | 1.19 | 1.19 | 1.19 | 1.19 |
| Set::Object | 1.42 | 1.42 | 1.42 | 1.42 | 1.42 |
| Set::Scalar | 1.29 | 1.29 | 1.29 | 1.29 | 1.29 |
| Signal::Mask | 0.008 | 0.008 | 0.008 | 0.008 | 0.008 |
| Smart::Comments | 1.06 | 1.06 | 1.06 | 1.06 | 1.06 |
| SOAP::Lite | 1.27 | 1.27 | 1.27 | 1.27 | 1.27 |
| SOAP::WSDL | n/a | n/a | 3.004 | n/a | n/a |
| Socket | n/a | 2.037 | n/a | n/a | n/a |
| Socket6 | 0.29 | 0.29 | 0.29 | 0.29 | 0.29 |
| Software::License | n/a | 0.104004 | 0.104002 | 0.104004 | 0.104004 |
| Software::License::CCpack | 1.11 | 1.11 | 1.11 | 1.11 | 1.11 |
| Sort::ArrayOfArrays | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| Sort::Key | 1.33 | 1.33 | 1.33 | 1.33 | 1.33 |
| Sort::Naturally | 1.03 | 1.03 | 1.03 | 1.03 | 1.03 |
| Sort::Versions | 1.62 | n/a | n/a | n/a | n/a |
| Specio | 0.47 | 0.48 | 0.48 | 0.48 | 0.48 |
| Spiffy | 0.46 | 0.46 | 0.46 | 0.46 | 0.46 |
| Spreadsheet::ParseExcel | 0.65 | 0.65 | 0.65 | 0.65 | 0.65 |
| Spreadsheet::ParseXLSX | 0.27 | 0.27 | 0.27 | 0.27 | 0.27 |
| Spreadsheet::WriteExcel | 2.40 | 2.40 | 2.40 | 2.40 | 2.40 |
| Spreadsheet::XLSX | 0.17 | 0.17 | 0.17 | 0.17 | 0.17 |
| Starlet | 0.31 | 0.31 | 0.31 | 0.31 | 0.31 |
| Statistics::Basic | 1.6611 | 1.6611 | 1.6611 | 1.6611 | 1.6611 |
| Statistics::Descriptive | 3.0800 | 3.0801 | 3.0800 | 3.0800 | 3.0801 |
| Statistics::Distributions | 1.02 | 1.02 | 1.02 | 1.02 | 1.02 |
| Statistics::Lite | 3.62 | n/a | n/a | n/a | n/a |
| Statistics::R | 0.34 | 0.34 | 0.34 | 0.34 | 0.34 |
| Statistics::TTest | 1.1 | 1.1 | 1.1 | 1.1 | 1.1 |
| Storable | n/a | 3.25 | n/a | n/a | n/a |
| Stream::Buffered | 0.03 | 0.03 | 0.03 | 0.03 | 0.03 |
| strictures | 2.000006 | 2.000006 | 2.000006 | 2.000006 | 2.000006 |
| String::Approx | 3.28 | 3.28 | 3.28 | 3.28 | 3.28 |
| String::HexConvert | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 |
| String::Print | 0.94 | 0.94 | 0.94 | 0.94 | 0.94 |
| String::Random | 0.32 | 0.32 | 0.32 | 0.32 | 0.32 |
| String::ShellQuote | n/a | 1.04 | 1.04 | 1.04 | 1.04 |
| Sub::Exporter | 0.989 | 0.990 | 0.988 | 0.989 | 0.990 |
| Sub::Exporter::ForMethods | 0.100054 | 0.100055 | 0.100054 | 0.100055 | 0.100055 |
| Sub::Exporter::Progressive | 0.001013 | 0.001013 | 0.001013 | 0.001013 | 0.001013 |
| Sub::HandlesVia | 0.016 | 0.050000 | 0.045 | 0.050000 | 0.050000 |
| Sub::Identify | 0.14 | 0.14 | 0.14 | 0.14 | 0.14 |
| Sub::Info | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 |
| Sub::Install | 0.929 | 0.929 | 0.928 | 0.929 | 0.929 |
| Sub::Name | 0.26 | 0.27 | 0.26 | 0.27 | 0.27 |
| Sub::Quote | 2.006006 | 2.006008 | 2.006006 | 2.006008 | 2.006008 |
| Sub::Uplevel | 0.2800 | 0.2800 | 0.2800 | 0.2800 | 0.2800 |
| SUPER | 1.20190531 | 1.20190531 | 1.20190531 | 1.20190531 | 1.20190531 |
| SVG | 2.86 | 2.87 | 2.87 | 2.87 | 2.87 |
| SVG::Graph | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 |
| SWF::File | 0.033 | 0.033 | 0.033 | 0.033 | 0.033 |
| Switch | 2.17 | 2.17 | 2.17 | 2.17 | 2.17 |
| Switch::Plain | 0.0501 | 0.0501 | 0.0501 | 0.0501 | 0.0501 |
| Symbol::Global::Name | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 |
| Sys::Info | 0.7811 | 0.7811 | 0.7811 | 0.7811 | 0.7811 |
| Sys::Info::Base | 0.7807 | 0.7807 | 0.7807 | 0.7807 | 0.7807 |
| Sys::Info::Driver::Linux | 0.7905 | 0.7905 | 0.7905 | 0.7905 | 0.7905 |
| Sys::RunAlone | 0.13 | 0.13 | 0.13 | 0.13 | 0.13 |
| Sys::SigAction | n/a | 0.23 | n/a | 0.23 | n/a |
| System::Timeout | 0.07 | 0.07 | 0.07 | 0.07 | 0.07 |
| Task::Weaken | 1.06 | 1.06 | 1.06 | 1.06 | 1.06 |
| Template | 3.010 | 3.101 | 3.101 | 3.101 | 3.101 |
| Template::Plugin::CGI | n/a | 3.101 | 3.101 | 3.101 | 3.101 |
| Term::ProgressBar | 2.22 | 2.23 | 2.23 | 2.23 | 2.23 |
| Term::ReadKey | 2.38 | 2.38 | 2.38 | 2.38 | 2.38 |
| Term::ReadLine | 1.17 | 1.12 | 1.17 | 1.17 | 1.17 |
| Term::Size::Any | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 |
| Term::Size::Perl | 0.031 | 0.031 | 0.031 | 0.031 | 0.031 |
| Term::Table | 0.016 | 0.016 | 0.016 | 0.016 | 0.016 |
| Term::UI | 0.50 | 0.50 | 0.50 | 0.50 | 0.50 |
| Test2::Plugin::NoWarnings | 0.09 | 0.09 | 0.09 | 0.09 | 0.09 |
| Test2::Suite | 0.000145 | 0.000155 | 0.000145 | 0.000155 | 0.000155 |
| Test2::Tools::Explain | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 |
| Test::API | n/a | 0.010 | n/a | 0.010 | 0.010 |
| Test::Base | 0.89 | 0.89 | 0.89 | 0.89 | 0.89 |
| Test::CheckDeps | 0.010 | 0.010 | 0.010 | 0.010 | 0.010 |
| Test::Class | 0.52 | 0.52 | 0.52 | 0.52 | 0.52 |
| Test::CleanNamespaces | 0.24 | 0.24 | 0.24 | 0.24 | 0.24 |
| Test::Compile | v3.0.1 | v3.3.1 | v3.1.0 | v3.2.2 | v3.3.1 |
| Test::CPAN::Meta | 0.25 | 0.25 | 0.25 | 0.25 | 0.25 |
| Test::CPAN::Meta::JSON | 0.16 | 0.16 | 0.16 | 0.16 | 0.16 |
| Test::Deep | 1.130 | 1.204 | 1.130 | 1.204 | 1.204 |
| Test::Differences | 0.69 | 0.70 | 0.69 | 0.69 | 0.70 |
| Test::EOL | 2.02 | 2.02 | 2.02 | 2.02 | 2.02 |
| Test::Exception | 0.43 | 0.43 | 0.43 | 0.43 | 0.43 |
| Test::Exports | n/a | n/a | n/a | 1 | 1 |
| Test::FailWarnings | 0.008 | 0.008 | 0.008 | 0.008 | 0.008 |
| Test::Fatal | 0.016 | 0.017 | 0.016 | 0.017 | 0.017 |
| Test::File | 1.992 | 1.993 | 1.992 | 1.993 | 1.993 |
| Test::File::ShareDir | 1.001002 | 1.001002 | 1.001002 | 1.001002 | 1.001002 |
| Test::Fork | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 |
| Test::Harness | n/a | 3.47 | n/a | n/a | n/a |
| Test::InDistDir | 1.112071 | 1.112071 | 1.112071 | 1.112071 | 1.112071 |
| Test::Inter | 1.09 | 1.10 | 1.09 | 1.10 | 1.10 |
| Test::JSON | 0.11 | 0.11 | 0.11 | 0.11 | 0.11 |
| Test::Kwalitee | 1.28 | 1.28 | 1.28 | 1.28 | 1.28 |
| Test::LeakTrace | 0.17 | 0.17 | 0.17 | 0.17 | 0.17 |
| Test::Lib | n/a | 0.003 | n/a | 0.003 | 0.003 |
| Test::LongString | 0.17 | 0.17 | 0.17 | 0.17 | 0.17 |
| Test::Memory::Cycle | 1.06 | 1.06 | 1.06 | 1.06 | 1.06 |
| Test::MockModule | 0.177.0 | 0.177.0 | 0.177.0 | 0.177.0 | 0.177.0 |
| Test::MockObject | 1.20200122 | 1.20200122 | 1.20200122 | 1.20200122 | 1.20200122 |
| Test::MockRandom | 1.01 | 1.01 | 1.01 | 1.01 | 1.01 |
| Test::MockTime | 0.17 | 0.17 | 0.17 | 0.17 | 0.17 |
| Test::MockTime::HiRes | n/a | 0.08 | n/a | 0.08 | 0.08 |
| Test::More::UTF8 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 |
| Test::Most | 0.37 | 0.38 | 0.38 | 0.38 | 0.38 |
| Test::Needs | 0.002009 | 0.002010 | 0.002009 | 0.002010 | 0.002010 |
| Test::NoTabs | 2.02 | 2.02 | 2.02 | 2.02 | 2.02 |
| Test::NoWarnings | 1.06 | 1.06 | 1.06 | 1.06 | 1.06 |
| Test::Number::Delta | 1.06 | 1.06 | 1.06 | 1.06 | 1.06 |
| Test::Object | 0.08 | 0.08 | 0.08 | 0.08 | 0.08 |
| Test::Output | 1.033 | 1.034 | 1.033 | 1.033 | 1.034 |
| Test::PAUSE::Permissions | 0.07 | 0.07 | 0.07 | 0.07 | 0.07 |
| Test::Pod | 1.52 | 1.52 | 1.52 | 1.52 | 1.52 |
| Test::Pod::Content | 0.0.6 | 0.0.6 | 0.0.6 | 0.0.6 | 0.0.6 |
| Test::Pod::Coverage | 1.10 | 1.10 | 1.10 | 1.10 | 1.10 |
| Test::Regexp | 2017040101 | 2017040101 | 2017040101 | 2017040101 | 2017040101 |
| Test::Reporter | 1.62 | 1.62 | 1.62 | 1.62 | 1.62 |
| Test::Reporter::Transport::Legacy | 1.59 | 1.59 | 1.59 | 1.59 | 1.59 |
| Test::Requires | 0.11 | 0.11 | 0.11 | 0.11 | 0.11 |
| Test::RequiresInternet | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 |
| Test::SharedFork | 0.35 | 0.35 | 0.35 | 0.35 | 0.35 |
| Test::Simple | 1.302190 | 1.302195 | 1.302191 | 1.302195 | 1.302195 |
| Test::SubCalls | 1.10 | 1.10 | 1.10 | 1.10 | 1.10 |
| Test::Sys::Info | 0.23 | 0.23 | 0.23 | 0.23 | 0.23 |
| Test::TCP | 2.22 | 2.22 | 2.22 | 2.22 | 2.22 |
| Test::Time | 0.08 | 0.092 | 0.092 | 0.092 | 0.092 |
| Test::Trap | 0.3.4 | 0.3.5 | 0.3.5 | 0.3.5 | 0.3.5 |
| Test::UseAllModules | 0.17 | 0.17 | 0.17 | 0.17 | 0.17 |
| Test::utf8 | 1.02 | 1.02 | 1.02 | 1.02 | 1.02 |
| Test::Warn | 0.36 | 0.37 | 0.37 | 0.37 | 0.37 |
| Test::Warnings | 0.031 | 0.031 | 0.031 | 0.031 | 0.031 |
| Test::Weaken | 3.022000 | 3.022000 | 3.022000 | 3.022000 | 3.022000 |
| Test::Without::Module | 0.20 | 0.21 | 0.21 | 0.21 | 0.21 |
| Test::XML | 0.08 | 0.08 | 0.08 | 0.08 | 0.08 |
| Test::YAML | 1.07 | 1.07 | 1.07 | 1.07 | 1.07 |
| Text::Aligner | 0.16 | 0.16 | 0.16 | 0.16 | 0.16 |
| Text::ASCIITable | 0.22 | 0.22 | 0.22 | 0.22 | 0.22 |
| Text::Autoformat | 1.75 | 1.75 | 1.75 | 1.75 | 1.75 |
| Text::CSV | 2.01 | 2.03 | 2.02 | 2.02 | 2.02 |
| Text::CSV\_XS | 1.47 | 1.51 | 1.48 | 1.50 | 1.50 |
| Text::Diff | n/a | 1.45 | 1.45 | 1.45 | 1.45 |
| Text::Format | 0.62 | 0.62 | 0.62 | 0.62 | 0.62 |
| Text::FormatTable | 1.03 | 1.03 | 1.03 | 1.03 | 1.03 |
| Text::Glob | n/a | 0.11 | 0.11 | 0.11 | 0.11 |
| Text::Haml | 0.990118 | 0.990118 | 0.990118 | 0.990118 | 0.990118 |
| Text::Levenshtein | 0.15 | 0.15 | 0.15 | 0.15 | 0.15 |
| Text::LevenshteinXS | 0.03 | 0.03 | 0.03 | 0.03 | 0.03 |
| Text::Password::Pronounceable | 0.30 | 0.30 | 0.30 | 0.30 | 0.30 |
| Text::Quoted | 2.10 | 2.10 | 2.10 | 2.10 | 2.10 |
| Text::Reform | 1.20 | 1.20 | 1.20 | 1.20 | 1.20 |
| Text::SimpleTable | 2.07 | 2.07 | 2.07 | 2.07 | 2.07 |
| Text::Soundex | 3.05 | 3.05 | 3.05 | 3.05 | 3.05 |
| Text::Table | 1.134 | 1.135 | 1.135 | 1.135 | 1.135 |
| Text::Template | n/a | 1.61 | 1.61 | 1.61 | 1.61 |
| Text::Template::Simple | 0.91 | 0.91 | 0.91 | 0.91 | 0.91 |
| Text::Unidecode | 1.30 | 1.30 | 1.30 | 1.30 | 1.30 |
| Text::vFile::asData | 0.08 | 0.08 | 0.08 | 0.08 | 0.08 |
| Text::WikiFormat | 0.81 | 0.81 | 0.81 | 0.81 | 0.81 |
| Text::Wrapper | 1.05 | 1.05 | 1.05 | 1.05 | 1.05 |
| Thread::Bless | 0.09 | 0.09 | 0.09 | 0.09 | 0.09 |
| Thread::Cancel | 1.13 | 1.13 | 1.13 | 1.13 | 1.13 |
| Thread::Cleanup | 0.07 | 0.07 | 0.07 | 0.07 | 0.07 |
| Thread::Conveyor | 0.20 | 0.20 | 0.20 | 0.20 | 0.20 |
| Thread::Conveyor::Monitored | 0.15 | 0.15 | 0.15 | 0.15 | 0.15 |
| Thread::Exit | 0.10 | 0.10 | 0.10 | 0.10 | 0.10 |
| Thread::Pool | 0.35 | 0.35 | 0.35 | 0.35 | 0.35 |
| Thread::Queue::Any | 1.16 | 1.16 | 1.16 | 1.16 | 1.16 |
| Thread::Rand | 0.08 | 0.08 | 0.08 | 0.08 | 0.08 |
| Thread::Running | 0.08 | 0.08 | 0.08 | 0.08 | 0.08 |
| Thread::Serialize | 1.02 | 1.02 | 1.02 | 1.02 | 1.02 |
| Thread::SigMask | 0.004 | 0.004 | 0.004 | 0.004 | 0.004 |
| threads::shared | n/a | 1.59 | n/a | n/a | n/a |
| Thread::Stack | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| Thread::State | 0.09 | 0.09 | 0.09 | 0.09 | 0.09 |
| Thread::Tie | 0.15 | 0.15 | 0.15 | 0.15 | 0.15 |
| Thread::Use | 0.06 | 0.06 | 0.06 | 0.06 | 0.06 |
| Throwable | 1.001 | 1.001 | n/a | 1.001 | 1.001 |
| Tie::IxHash | 1.23 | 1.23 | 1.23 | 1.23 | 1.23 |
| Tie::ToObject | 0.03 | 0.03 | 0.03 | 0.03 | 0.03 |
| Time::Duration | 1.21 | 1.21 | 1.21 | 1.21 | 1.21 |
| Time::Duration::Parse | 0.16 | 0.16 | 0.16 | 0.16 | 0.16 |
| Time::HiRes | 1.9764 | 1.9764 | n/a | n/a | n/a |
| Time::Local | n/a | 1.35 | n/a | n/a | n/a |
| Time::NT | 0.007 | 0.007 | 0.007 | 0.007 | 0.007 |
| Time::ParseDate | 2015.103 | 2015.103 | 2015.103 | 2015.103 | 2015.103 |
| Time::Piece | 1.3401 | n/a | n/a | n/a | n/a |
| Tk | 804.036 | 804.036 | 804.036 | 804.036 | 804.036 |
| Tree::DAG\_Node | 1.32 | 1.32 | 1.32 | 1.32 | 1.32 |
| Tree::Simple | 1.34 | 1.34 | 1.34 | 1.34 | 1.34 |
| true | n/a | n/a | n/a | v1.0.2 | v1.0.2 |
| Try::Tiny | 0.31 | 0.31 | 0.31 | 0.31 | 0.31 |
| Type::API | n/a | 1.001 | 1.001 | 1.001 | 1.001 |
| Types::Path::Tiny | 0.006 | 0.006 | 0.006 | 0.006 | 0.006 |
| Types::Serialiser | 1.01 | 1.01 | 1.01 | 1.01 | 1.01 |
| Type::Tie | 0.015 | n/a | n/a | n/a | n/a |
| Type::Tiny | 1.012004 | 2.004000 | 1.016010 | 2.004000 | 2.004000 |
| Type::Tiny::XS | 0.022 | 0.025 | 0.025 | 0.025 | 0.025 |
| Unicode::Collate | n/a | 1.31 | n/a | n/a | n/a |
| Unicode::LineBreak | 2019.001 | 2019.001 | 2019.001 | 2019.001 | 2019.001 |
| Unicode::UTF8 | 0.62 | 0.62 | 0.62 | 0.62 | 0.62 |
| UNIVERSAL::can | 1.20140328 | 1.20140328 | 1.20140328 | 1.20140328 | 1.20140328 |
| UNIVERSAL::isa | 1.20171012 | 1.20171012 | 1.20171012 | 1.20171012 | 1.20171012 |
| UNIVERSAL::ref | 0.14 | 0.14 | 0.14 | 0.14 | 0.14 |
| UNIVERSAL::require | 0.19 | 0.19 | 0.19 | 0.19 | 0.19 |
| Unix::Processors | 2.046 | 2.046 | 2.046 | 2.046 | 2.046 |
| URI | 5.10 | 5.21 | 5.17 | 5.19 | 5.19 |
| URI::Escape::XS | 0.14 | 0.14 | 0.14 | 0.14 | 0.14 |
| Variable::Magic | 0.62 | 0.63 | 0.63 | 0.63 | 0.63 |
| version | 0.9929 | 0.9929 | 0.9929 | 0.9929 | 0.9929 |
| Want | 0.29 | 0.29 | 0.29 | 0.29 | 0.29 |
| WWW::Curl | 4.17 | 4.17 | 4.17 | 4.17 | 4.17 |
| WWW::Form::UrlEncoded | 0.26 | 0.26 | 0.26 | 0.26 | 0.26 |
| WWW::RobotRules | 6.02 | 6.02 | 6.02 | 6.02 | 6.02 |
| XML::CommonNS | 0.06 | 0.06 | 0.06 | 0.06 | 0.06 |
| XML::Compile | 1.63 | 1.63 | 1.63 | 1.63 | 1.63 |
| XML::Compile::C14N | 0.95 | 0.95 | 0.95 | 0.95 | 0.95 |
| XML::Compile::Cache | 1.06 | 1.06 | 1.06 | 1.06 | 1.06 |
| XML::Compile::Dumper | 0.14 | 0.14 | 0.14 | 0.14 | 0.14 |
| XML::Compile::SOAP | 3.27 | 3.28 | 3.28 | 3.28 | 3.28 |
| XML::Compile::SOAP::Daemon | 3.14 | 3.15 | 3.14 | 3.15 | 3.15 |
| XML::Compile::Tester | 0.91 | 0.91 | 0.91 | 0.91 | 0.91 |
| XML::Compile::WSA | 0.95 | 0.95 | 0.95 | 0.95 | 0.95 |
| XML::Compile::WSDL11 | 3.08 | 3.08 | 3.08 | 3.08 | 3.08 |
| XML::Compile::WSS | 1.14 | 1.14 | 1.14 | 1.14 | 1.14 |
| XML::Compile::WSS::Signature | 2.02 | 2.02 | 2.02 | 2.02 | 2.02 |
| XML::DOM | 1.46 | 1.46 | 1.46 | 1.46 | 1.46 |
| XML::DOM::XPath | 0.14 | 0.14 | 0.14 | 0.14 | 0.14 |
| XML::Filter::BufferText | 1.01 | 1.01 | 1.01 | 1.01 | 1.01 |
| XML::LibXML | 2.0207 | 2.0209 | 2.0208 | 2.0208 | 2.0209 |
| XML::LibXML::Simple | 1.01 | 1.01 | 1.01 | 1.01 | 1.01 |
| XML::Namespace | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 |
| XML::NamespaceFactory | 1.02 | 1.02 | 1.02 | 1.02 | 1.02 |
| XML::NamespaceSupport | 1.12 | 1.12 | 1.12 | 1.12 | 1.12 |
| XML::Parser | 2.36 | 2.36 | 2.36 | 2.36 | 2.36 |
| XML::Parser::Lite | 0.722 | 0.722 | 0.722 | 0.722 | 0.722 |
| XML::RegExp | 0.04 | 0.04 | undef | 0.04 | 0.04 |
| XML::RSS | 1.62 | 1.62 | 1.62 | 1.62 | 1.62 |
| XML::SAX | 1.02 | 1.02 | 1.02 | 1.02 | 1.02 |
| XML::SAX::Base | 1.09 | 1.09 | 1.09 | 1.09 | 1.09 |
| XML::SAX::Expat | 0.51 | 0.51 | 0.51 | 0.51 | 0.51 |
| XML::SAX::Writer | 0.57 | 0.57 | 0.57 | 0.57 | 0.57 |
| XML::SemanticDiff | 1.0007 | 1.0007 | 1.0007 | 1.0007 | 1.0007 |
| XML::Simple | 2.25 | 2.25 | 2.25 | 2.25 | 2.25 |
| XML::Twig | 3.52 | 3.52 | 3.52 | 3.52 | 3.52 |
| XML::Writer | 0.900 | 0.900 | 0.900 | 0.900 | 0.900 |
| XML::XML2JSON | 0.06 | 0.06 | 0.06 | 0.06 | 0.06 |
| XML::XPath | 1.44 | 1.48 | 1.48 | 1.48 | 1.48 |
| XML::XPathEngine | 0.14 | 0.14 | 0.14 | 0.14 | 0.14 |
| XSLoader | 0.24 | 0.24 | n/a | n/a | n/a |
| XS::Parse::Keyword | n/a | 0.38 | n/a | 0.33 | 0.36 |
| XS::Parse::Sublike | n/a | 0.18 | n/a | 0.17 | 0.18 |
| XString | 0.005 | 0.005 | 0.005 | 0.005 | 0.005 |
| XXX | 0.38 | 0.38 | 0.38 | 0.38 | 0.38 |
| YAML | n/a | 1.30 | 1.30 | 1.30 | 1.30 |
| YAML::LibYAML | 0.83 | 0.88 | 0.85 | 0.88 | 0.88 |
| YAML::PP | 0.032 | 0.036 | 0.035 | 0.036 | 0.036 |
| YAML::Syck | 1.34 | 1.34 | 1.34 | 1.34 | 1.34 |
| YAML::Tiny | 1.73 | 1.74 | 1.73 | 1.74 | 1.74 |


 

Personal installation using cpanm
[top](#top)
Perl modules can be installed locally in your /home directory easily using
[cpanm](https://metacpan.org/dist/App-cpanminus/view/bin/cpanm).


* First, allocate an interactive session with plenty of memory:
 
```
*user*@biowulf:~$ sinteractive --mem=20g
```
* Choose a location. Most likely ~/perl will work just fine.
* Append the following lines to your ~/.bashrc file, using the following heredoc:
 
```
*user*@*node*:~$ cat << EOF >> ~/.bashrc

#------------------------------------------------------------------------------
# For local perl
#------------------------------------------------------------------------------
export LOCALPERL=~/perl
export PERL5LIB=$LOCALPERL:$LOCALPERL/lib/perl5
export PERL_CPANM_HOME=$LOCALPERL/cpanm
export PERL_CPANM_OPT="-l $LOCALPERL"
export PATH=$LOCALPERL/bin:$PATH
export PERLDOC_PAGER="less -isR"
EOF
```
* Source your ~/.bashrc file:
 
```
*user*@biowulf:~$ source ~/.bashrc
```
* Create local directories:
 
```
*user*@biowulf:~$ mkdir -p $LOCALPERL
*user*@biowulf:~$ mkdir -p $PERL_CPANM_HOME
```
* Install a module:
 
```
*user*@biowulf:~$ cpanm TableData
```
* Read about the module:
 
```
*user*@biowulf:~$ perldoc TableData
```
* See where it lives:
 
```
*user*@biowulf:~$ perldoc -l TableData
~/perl/lib/perl5/TableData.pm
```
* More detail:
 
```
*user*@biowulf:~$ tree -L 3 ~/perl
~/perl
├── bin
├── cpanm
│   ├── build.log -> ~/perl/cpanm/work/1686253267.119612/build.log
│   ├── latest-build -> ~/perl/cpanm/work/1686253267.119612
│   └── work
│       └── 1686253267.119612
├── lib
│   └── perl5
│       ├── 5.26.3
│       ├── Role
│       ├── RoleBundle
│       ├── TableData
│       ├── TableData.pm
│       ├── TableDataRole
│       └── x86_64-linux-thread-multi
└── man
    └── man3
        ├── RoleBundle::TinyCommons::Iterator.3pm
        ├── Role::TinyCommons::Iterator::Basic.3pm
        ├── Role::TinyCommons::Iterator::Bidirectional.3pm
        ├── Role::TinyCommons::Iterator::Circular.3pm
        ├── Role::TinyCommons::Iterator::Resettable.3pm
        ├── TableData.3pm
        ├── TableDataRole::Spec::Basic.3pm
        └── TableData::Test::Spec::Basic.3pm
```


Some modules may require specific libraries or executables. Other
modules may have trouble compiling. If you have difficulties, please send an
email to staff@hpc.nih.gov.




