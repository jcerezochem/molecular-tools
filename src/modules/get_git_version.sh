#!/bin/bash

git_hash=$(git describe --long --dirty --always)
git_date=$(git show -s --format=%ci)
user=$USER

cat <<EOF > version.f90
    subroutine print_version()

        write(0,'(/A)') "GIT VERSION INFO"
        write(0,'(A)')   " Commit id  : $git_hash"
        write(0,'(A)')   " Commit date: $git_date"
        
        return
    end subroutine 
EOF
