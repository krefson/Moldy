
#
# Extension to LaTeX2HTML supply support for the "superfig"
# LaTeX style, used in Moldy.
# Change Log:
# ===========
 
package main;
#
#  Make the floatingfigure environment be translated as
#  an ordinary figure, ignoring the mandatory width.
#
#
 
sub do_env_partfigure {
    local($_) = @_;
 
    $contents =~ s/$optional_arg_rx//o;    # ditch [tbp]
    &process_environment("figure", $global{'max_id'}++);
    }

sub do_env_Argdescription {
    local($_) = @_;
    #RRM - catch nested lists
    $_ = &translate_environments($_);
    $* = 1;                     # Multiline matching ON
    s/$item_description_rx/<DT><TT>$1<\/TT>\n<DD>/g;
    $* = 0;                     # Multiline matching OFF
    &do_env_description($_, " COMPACT");
}

sub do_env_Litdescription {
    local($_) = @_;
    #RRM - catch nested lists
    $_ = &translate_environments($_);
    $* = 1;                     # Multiline matching ON
    s/$item_description_rx/<DT><TT>$1<\/TT>\n<DD>/g;
    $* = 0;                     # Multiline matching OFF
    &do_env_description($_, " COMPACT");
}

sub do_env_Fndescription {
    local($_) = @_;
    #RRM - catch nested lists
    $_ = &translate_environments($_);
    $* = 1;                     # Multiline matching ON
    s/$item_description_rx/<DT><I>$1<\/I>\n<DD>/g;
    $* = 0;                     # Multiline matching OFF
    &do_env_description($_, " COMPACT");
}

sub do_cmd_Emph { &styled_text_chunk('B','em','font','variant','','', @_); }
sub do_cmd_Lit { &styled_text_chunk('TT','em','font','variant','','', @_); }
sub do_cmd_Fname { &styled_text_chunk('I','em','font','variant','','', @_); }

&process_commands_in_tex (<<_RAW_ARG_CMDS_);
partfigure # []
_RAW_ARG_CMDS_
 
&ignore_commands( <<_IGNORED_CMDS_);
superfig
saferagged
_IGNORED_CMDS_
1;                              # This must be the last line
