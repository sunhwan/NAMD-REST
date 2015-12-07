
tclForcesScript {
    addatom 1
    addatom 2
    addatom 3
    addatom 4
    addatom 5
    addatom 6

proc calcforces {} {
set ts 0

    #loadcoords xyz
    #foreach {atomid coord} [array get xyz] {
    #  puts "XYZ:$atomid $ts $coord"
    #}

    loadforces f
    foreach {atomid force} [array get f] {
       puts "FORCEext:$atomid $ts  $force"
    }

    loadtotalforces ft
    foreach {atomid force} [array get ft] {
       puts "FORCEtot:$atomid $ts  $force"
    }

    print 'xyz'
}
}

