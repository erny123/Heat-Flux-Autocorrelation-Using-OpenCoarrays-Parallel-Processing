program coarraytest
    implicit none
    real:: a,b,c,d, sum,var,mean, dp[*]
    integer:: i,j, jsize,beginning, rate, end, end1, proc_num
    real, dimension(20000000):: Jx, Jxm !, corr
    real, dimension(1000000):: corr
    integer:: skip_lines = 4
    call system_clock(beginning, rate)
    !get # of processes
    proc_num = num_images()
    print*, this_image(), ": ", this_image()
    !begin reading data
    open(10, file='DiamHeat.log', status='old')
    do i = 1,skip_lines
        read(10,*)
    end do
    do i = 1, 20000000
        read(10,*) a, b, Jx(i), c, d
    end do
    call system_clock(end)
    print *, "elapsed time for reading: ", real(end - beginning) / real(rate)
    close(10)
    !finished reading data
    jsize = size(Jx)

    !calculate mean
    mean = sum(Jx)/jsize
    Jxm(:) = Jx(:)-mean

    !calculate variance
    var = dot_product(Jxm,Jxm)/jsize

    !open to write file
    open(20, file='acorr_coarray2.txt')
    call system_clock(end1)
    !print *, "Size of Jx and time: "
    !print *, jsize, real(end1 - beginning) / real(rate)
    sync all
    
    !print *, "this is var: ", var

    do j = 0, 2000000/proc_num
    
        dp = dot_product(Jxm(this_image()+proc_num*j:jsize),Jxm(1:jsize+1-this_image()-proc_num*j))/var &
        /(jsize+1-this_image()-proc_num*j)
        sync all       
        !print *, this_image(), ": ", (jsize+1-this_image()-proc_num*j)
        if(this_image() == 1) then
            !do i = 1, proc_num
                !acorr(i+j*proc_num) = dp[i]
                !corr(i+j*proc_num) = dp[i]
		!write(20,*) dp[
            !end do
	    write(20,*) (dp[i], i = 1,proc_num)
            if(j == 1) then
                call system_clock(end)
                print *, "correlation time magnitude 1e0 elapsed time: ", real(end - beginning) / real(rate)
            else if(j == 10) then
                call system_clock(end)
                print *, "correlation time magnitude 1e1 elapsed time: ", real(end - beginning) / real(rate)
            else if(j == 100) then
                call system_clock(end)
                print *, "correlation time magnitude 1e2 elapsed time: ", real(end - beginning) / real(rate)
            else if(j == 1000) then
                call system_clock(end)
                print *, "correlation time magnitude 1e3 elapsed time: ", real(end - beginning) / real(rate)
            else if(j == 10000) then
                call system_clock(end)
                print *, "correlation time magnitude 1e4 elapsed time: ", real(end - beginning) / real(rate)
            else if(j == 100000) then
                call system_clock(end)
                print *, "correlation time magnitude 1e5 elapsed time: ", real(end - beginning) / real(rate)
            else if(j == 1000000) then
                call system_clock(end)
                print *, "correlation time magnitude 1e6 elapsed time: ", real(end - beginning) / real(rate)
            else if(j == 10000000) then
                call system_clock(end)
                print *, "correlation time magnitude 1e7 elapsed time: ", real(end - beginning) / real(rate)
            end if 

            !print*, acorr
        end if 


    end do 
    sync all
    !if(this_image()==1) then 
        !print*, size(acorr)
        !print*, acorr
        !write(20,*) corr
    !end if 
call system_clock(end)
print *, "total time: ", real(end - beginning) / real(rate)


    close(20)

end program
