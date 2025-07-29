function p =polymin(p1,p2) 

p = [zeros(1, size(p2,2)-size(p1,2)) p1]-[zeros(1, size(p1,2)-size(p2,2)) p2] ;