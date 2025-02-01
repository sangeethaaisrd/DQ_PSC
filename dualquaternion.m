classdef dualquaternion
  
  properties
    qr (1,1) quaternion;  % Real quaternion (4x1 vector)
    qd (1,1) quaternion;  % dual quaternion (4x1 vector)
  end

  methods

    function dq = dualquaternion(qr,qd)
      % Constructor
      if nargin == 0
        dq.qr = quaternion(0,0); % initializes real quaternion
        dq.qd = quaternion(0,0); % initializes dual quaternion
      elseif nargin == 1
          if isscalar(qr)
              dq.qr = quaternion(qr);
              dq.qd = quaternion();
          else
              dq.qr = quaternion(qr(1),qr(2:4));
              dq.qd = quaternion(qr(5),qr(6:8));
          end
      else
        dq.qr = qr;
        dq.qd = qd;
      end
    end

    function dq_out = uplus(dq)
      % Dual Quaternion unary addition
      dq_out = dq;
    end

    function dq_out = plus(dq1, dq2)
      % Quaternion addition
      dq_out = dualquaternion();
      if isa(dq1,'double')
          dq1 = dualquaternion(quaternion(dq1,[0 0 0]),quaternion());
      elseif isa(dq2,'double')
          dq2 = dualquaternion(quaternion(dq2,[0 0 0]),quaternion());
      end
      dq_out.qr = dq1.qr + dq2.qr;
      dq_out.qd = dq1.qd + dq2.qd;
    end

    function dq_out = uminus(dq)
      % Dual Quaternion unary substraction
      dq_out = dualquaternion();
      dq_out.qr = -dq.qr;
      dq_out.qd = -dq.qd;
    end

    function dq_out = minus(dq1,dq2)
      % Quaternion subtraction
      dq_out = dualquaternion();
      if isa(dq1,'double')
          dq1 = dualquaternion(quaternion(dq1,[0 0 0]),quaternion());
      elseif isa(dq2,'double')
          dq2 = dualquaternion(quaternion(dq2,[0 0 0]),quaternion());
      end
      dq_out.qr = dq1.qr - dq2.qr;
      dq_out.qd = dq1.qd - dq2.qd;
    end    
            
    function dq_out = mtimes(dq1, dq2)
       % Dual quaternion multiplication
       dq_out = dualquaternion();
       if isscalar(dq1)&& isscalar(dq2)
        if isa(dq1,'double')
              dq1 = dualquaternion(quaternion(dq1,[0 0 0]),quaternion());
        elseif isa(dq2,'double')
              dq2 = dualquaternion(quaternion(dq2,[0 0 0]),quaternion());
        end
        dq_out.qr = dq1.qr*dq2.qr;  
        dq_out.qd = dq1.qr*dq2.qd + dq1.qd*dq2.qr;
       else
           if ~isscalar(dq1)
                dq_out = dualquaternion(dq1*dq2vec(dq2));
           else
                error('Right matrix multiplication not defined!');
           end
     end
    end

    function dq_out = cross(dq1, dq2)
      % Dual quaternion multiplication
      dq_out = dualquaternion();
      dq_out.qr =  cross(dq1.qr,dq2.qr);
      dq_out.qd =  cross(dq1.qr,dq2.qd) + cross(dq1.qd,dq2.qr);
     end   

     function dq_out = conj(dq)
         % conjugate of a dual quaternion
      dq_out = dualquaternion();
      dq_out.qr =  conj(dq.qr);
      dq_out.qd =  conj(dq.qd);
     end            

     function dq_out = norm(dq)
         % norm of dual quaternion is a dual number; representing it as a
         % dual quaternion
         dq_out = dualquaternion();
         dq_out.qr = quaternion(norm(dq.qr));
         dq_out.qd = quaternion(dot(dq.qr,dq.qd)/norm(dq.qr));
     end

     function dq_out = inv(dq)
      % Dual quaternion inverse
      dq_out = dualquaternion();
      dq_out.qr = inv(dq.qr);
      dq_out.qd = -inv(dq.qr)*dq.qd*inv(dq.qr);
     end

     function dq_out = normalize(dq)
         % normalize dual quaternion
         dq_out = dq/norm(dq);
     end

     function dq_out = mrdivide(dq1,dq2)
      % Quaternion division
      if isa(dq1,'double')
          dq1 = dualquaternion(quaternion(dq1,[0 0 0]),quaternion());
      elseif isa(dq2,'double')
          dq2 = dualquaternion(quaternion(dq2,[0 0 0]),quaternion());
      end
      dq_out = dq1*inv(dq2);
     end    
    
    function dq_vec = dq2vec(dq)
        dq_vec = [dq.qr.s; dq.qr.v; dq.qd.s; dq.qd.v];
    end

    function r_b= rb_from_dq(dq)
        r_b = 2*conj(dq.qr)*dq.qd;
    end

    function r_a= ra_from_dq(dq)
        r_a = 2*dq.qd*conj(dq.qr);
    end   

    function [theta,d,I,m] = con_2_screw(dq1)  % converting DQ to screw parameters
       
        theta = 2*acos(dq1.qr.s);
        d = -2 * dq1.qd.s *(1/sqrt(dot(dq1.qr.v,dq1.qr.v)));
        I = dq1.qr.v*(1/sqrt(dot(dq1.qr.v,dq1.qr.v)));
        m = (dq1.qd.v - 0.5*I*d*dq1.qr.s)*(1/sqrt(dot(dq1.qr.v,dq1.qr.v)));
    end 


    function pow_dq = pow_fcn(dq1,n) % converting to DQ from screw . n is the power.
        % for n =1, direct conversion from screw to DQ
      [theta,d,I,m] =  con_2_screw(dq1); 
      s1 = cos(0.5*n*theta);
      v1 = I*sin(0.5*n*theta);
      s2 = -0.5*n*d*sin(0.5*n*theta);
      v2 = sin(0.5*n*theta)*m+ 0.5*n*d*cos(0.5*n*theta)*I;
      qr =quaternion(s1,v1);
      qd = quaternion(s2,v2);
      pow_dq = dualquaternion(qr,qd);

    end 

    function interp_dq= dqinterp(t,dq1,dq2)
        interp_mag = conj(dq1)*(dq2);
        interp_dq = dq1*pow_fcn(interp_mag,t);
    end 


    function dq_pdt_elementwise= times(dq1,dq2)
         s1 = size(dq1);
         s2 = size(dq2);
        
        if s1(2) == s2(2)

                for i = 1:1:s1(2)                    
                    dq_pdt_elementwise(1,i) = dq1(1,i)*dq2(1,i);
                end                
        elseif s1(1)==8 
                 for i = 1:1:s2(2)
                     dq_pdt_elementwise(1,i) = dq1*dq2(1,i);
                 end

        else 

                error('dimensions must match!');
        end
    end 


    function dq_cross_elementwise= cross_times(dq1,dq2)
         s1 = size(dq1);
         s2 = size(dq2);
        
        if s1(2) == s2(2)

                for i = 1:1:s1(2)                    
                    dq_cross_elementwise(1,i) = dq1(1,i)*dq2(1,i);
                end                
        else 

                error('dimensions must match!');
        end

   
    end


    function dq_out = minus_times(dq1,dq2)
      % Quaternion subtraction
        
        s1 = size(dq1);
        s2 = size(dq2);

        if size(s1) == size(s2)

                for i = 1:1:s1(2) 
                    dq_out(1,i) = dq1(1,i)-dq2(1,i);
                end                
        else 

                error('dimensions must match!');
        end
    end 
  end
end 



