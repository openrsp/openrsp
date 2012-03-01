!
! WAVSPH Solvent sphere radii for nuclei
! WAVDEN Density matrix to be used in potential evaluations
! WAVFCK Solvation contribution to the Fock matrix
! WAVOE  Orbital energies
! WAVPCM .TRUE. if wavelet-pcm is to be used
!
      LOGICAL WAVPCM, WAVPCN
      COMMON /WAVPC/ WAVSPH(MXCENT),                                    &
     &                WAVDEN(MXORBT), WAVFCK(MXORBT),                   &
     &                WAVOE(MXCORB), WAVSCFE,                           &
     &                WAVPCM, WAVPCN

