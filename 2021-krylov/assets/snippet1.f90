znorm = dnrm2( s,xr,1 )
vm(:,1) = xr / znorm
! Comienzan las iteraciones de Lanczos
do i = 1,m
    call dgemv( 'n',s,s,1.0_dp,sigma,s,vm(:,i),1,0.0_dp,w,1 )
    if ( i > 1 ) then
        w = w - ( h(i-1,i) * vm(:,i-1) )
    end if
    k = size(w, 1)
    h(i,i) = ddot( k,w,1,vm(:,i),1 )
    if ( i < m ) then
        w = w - ( h(i,i) * vm(:,i) )
        k = size(w, 1)
        h(i+1,i) = dnrm2( k,w,1 )
        h(i,i+1) = h(i+1,i)
        vm(:,i+1) = w / h(i+1,i)
    end if
end do
call sqrt_matrix( h,temp )
call dgemv( 'N',m,m,1.0_dp,temp,m,eid,1,0.0_dp,v,1 )
call dgemv( 'N',s,m,1.0_dp,vm,s,v,1,0.0_dp,res,1 )
r = znorm * res
