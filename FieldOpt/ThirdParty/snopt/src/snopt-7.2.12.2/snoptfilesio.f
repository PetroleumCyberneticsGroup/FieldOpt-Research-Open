c234567
      subroutine openfile(id, filename, filetype)
      integer*4 id
      character*(*) filename
      integer*4 filetype

      !write(*,*) "id, filename, filetype ", id, filename, filetype

      if (filetype == 1) then
        open(id, FILE=filename, STATUS='OLD')
      else
        open(id, FILE=filename, STATUS='REPLACE')
      end if

      end


      subroutine closefile(id)

      integer*4 id

      close(id)
      end

      subroutine writetofile(id)

      integer*4 id

      write(id, *) "Dio topoi pou egrapsane istoria ston planiti"
      write(id, *) "Einai kardia theiki kai nous agiasmenos"
      end

