      parameter (Lpowmax=6)
      dimension ixyzpow(3*(Lpowmax+1)*(Lpowmax+1))
cbs   the ones and zeros stand four odd and even powers of x,y,z
cbs   if you want to go higher than l=6, you have to look up
cbs   the powers yourself, and add them to the table
      data ixyzpow /
     *0,0,0,                                 ! s-function
     *0,1,0, 0,0,1, 1,0,0,                   ! p-functions
     *1,1,0, 0,1,1, 0,0,0,  1,0,1, 0,0,0,    ! d-functions
     *0,1,0, 1,1,1, 0,1,0,  0,0,1, 1,0,0,    ! f-functions
     *0,0,1, 1,0,0,                          ! f-functions
     *1,1,0, 0,1,1, 1,1,0,  0,1,1, 0,0,0,    ! g-functions
     *1,0,1, 0,0,0, 1,0,1,  0,0,0,           ! g-functions
     *0,1,0, 1,1,1, 0,1,0,  1,1,1, 0,1,0,    ! h-functions
     *0,0,1, 1,0,0, 0,0,1,  1,0,0, 0,0,1,    ! h-functions
     *1,0,0,                                 ! h-functions
     *1,1,0, 0,1,1, 1,1,0, 0,1,1, 1,1,0,     ! i-functions
     *0,1,1, 0,0,0, 1,0,1, 0,0,0, 1,0,1,     ! i-functions
     *0,0,0, 1,0,1, 0,0,0                    ! i-functions
     */
