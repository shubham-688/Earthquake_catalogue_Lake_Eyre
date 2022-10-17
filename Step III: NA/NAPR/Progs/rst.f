        character*8 cname
c
        write(6,*) 'model name'
        read(5,*) cname
        write(6,*) cname
c
        call remodl(cname)
c
        call setbrn(cname)
c
c       call ttimes
        call tqresidm(cname)
c
        stop
        end
