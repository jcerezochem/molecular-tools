#!/bin/bash

git_hash=$(git describe --long --dirty --always)
git_date=$(git show -s --format=%ci)
user=$USER

if [ -f version.f90 ]; then
    old_git_hash=$(grep "Commit id  :" version.f90)
    old_git_hash=${old_git_hash##*Commit id  : }
    old_git_hash=${old_git_hash%\"}
else 
    old_git_hash=""
fi


if [ "$old_git_hash" != "$git_hash" ]; then
    echo ""
    echo "*************************************************"
    echo "Git hash changed from previous compilation: "
    echo " Old: $old_git_hash"
    echo " New: $git_hash"
    echo " A new version.f90 file will be generated"
    echo "*************************************************"
    echo ""
cat <<EOF > version.f90
    subroutine print_version()

        use io, only: uout

        write(uout,'(/A)') "GIT VERSION INFO"
        write(uout,'(A)')   " Commit id  : $git_hash"
        write(uout,'(A)')   " Commit date: $git_date"
        
        return
    end subroutine 
EOF

else 
    echo ""
    echo "*************************************************"
    echo "Git hash did not change from previous compilation"
    echo "*************************************************"
    echo ""
fi
