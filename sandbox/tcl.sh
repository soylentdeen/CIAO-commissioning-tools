#! /usr/bin/env tclsh

foreach d [concat [
    list $tcl_library [lindex $tcl_pkgPath 0]] $auto_path [
         list [file dirname $tcl_library] [
         file dirname [lindex $tcl_pkgPath 0]] [
         file dirname [file dirname $tcl_library]] [
         file dirname [file dirname [lindex $tcl_pkgPath 0]]] [
         file dirname [file dirname [file dirname $tcl_library]]] [
         file dirname [file dirname [
              file dirname [lindex $tcl_pkgPath 0]]]]]] {

     if {[file exists $d/tclConfig.sh]} {
         puts $d/tclConfig.sh
         exit
     }
}

puts none
