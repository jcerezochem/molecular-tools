#!/bin/bash

git_hash=$(git describe --long --dirty --always)
git_date=$(git show -s --format=%ci)
user=$USER

cat <<EOF > version.f90
    subroutine print_version()

        write(6,'(/A)') "GIT VERSION INFO"
        write(6,'(A)')   " Commit id  : $git_hash"
        write(6,'(A)')   " Commit date: $git_date"
        
        return
    end subroutine 
EOF
