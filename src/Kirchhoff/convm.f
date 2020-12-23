      SUBROUTINE CONV( H, L, X, M, Y, N )                                       
      DIMENSION H(*), X(*), Y(*)                                                
      N=L+M-1                                                                   
      IF ( L.GT.M) GO TO 300                                                    
      YY=0.0                                                                    
      DO 50 I=1, N                                                              
      IF ( I. GT. L ) GO TO 301                                                 
      YY=0.0                                                                    
      DO 51 J=1, I                                                              
      K=I-J+1                                                                   
51    YY=YY+H(J)*X(K)                                                           
      Y(I)=YY                                                                   
      GO TO 50                                                                  
301   CONTINUE                                                                  
      IF( I.GT. M ) GO TO 310                                                   
      YY=0.0                                                                    
      DO 52 J=1, L                                                              
      K=I-J+1                                                                   
52    YY=YY+H(J)*X(K)                                                           
      Y(I)=YY                                                                   
      GO TO 50                                                                  
310   CONTINUE                                                                  
      YY=0.0                                                                    
      J1=I-M+1                                                                  
      DO 60 J=J1, L                                                             
      K=I-J+1                                                                   
      YY=YY+H(J)*X(K)                                                           
60    CONTINUE                                                                  
      Y(I)=YY                                                                   
50    CONTINUE                                                                  
      RETURN                                                                    
300   CONTINUE                                                                  
      DO 53  I=1, N                                                             
      IF ( I.GT. M) GO TO 302                                                   
      YY=0.0                                                                    
      DO 54 J=1, I                                                              
      K=I-J+1                                                                   
54    YY=YY+H(J)*X(K)                                                           
      Y(I)=YY                                                                   
      GO TO 53                                                                  
302   CONTINUE                                                                  
      IF ( I.GT. L ) GO TO 312                                                  
      YY=0.0                                                                    
      J1=I-M+1                                                                  
      DO 55 J=J1, I                                                             
      K=I-J+1                                                                   
55    YY=YY+H(J)*X(K)                                                           
      Y(I)=YY                                                                   
      GO TO 53                                                                  
312   CONTINUE                                                                  
      YY=0.0                                                                    
      J1=I-M+1                                                                  
      DO 65 J=J1,L                                                              
      K=I-J+1                                                                   
      YY=YY+H(J)*X(K)                                                           
65    CONTINUE                                                                  
      Y(I)=YY                                                                   
53    CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
