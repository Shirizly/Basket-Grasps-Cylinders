function Avector = EEheightAvector(dp1,dp2,nj1,nj2,nk1,nk2,vj1,vj2,vk1,vk2)
%EEHEIGHTAVECTOR
%    AVECTOR = EEHEIGHTAVECTOR(DP1,DP2,NJ1,NJ2,NK1,NK2,VJ1,VJ2,VK1,VK2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    02-Apr-2021 01:49:19

t2 = dp1.*nj1;
t3 = dp1.*nj2;
t4 = dp2.*nj1;
t5 = dp2.*nj2;
t6 = dp1.*nk1;
t7 = dp1.*nk2;
t8 = dp2.*nk1;
t9 = dp2.*nk2;
t10 = nj1.*nk2;
t11 = nj2.*nk1;
t12 = nj1.*vj1;
t13 = nj2.*vj2;
t14 = nj1.*vk1;
t15 = nk1.*vj1;
t16 = nj2.*vk2;
t17 = nk2.*vj2;
t18 = nk1.*vk1;
t19 = nk2.*vk2;
t20 = -t4;
t21 = -t8;
t22 = -t11;
t23 = -t14;
t24 = -t16;
t25 = -t18;
t26 = -t19;
t27 = t2+t5;
t28 = t6+t9;
t29 = t3+t20;
t30 = t7+t21;
t31 = t10+t22;
t33 = t12+t13+t23+t24;
t34 = t15+t17+t25+t26;
t32 = 1.0./t31;
t35 = (nk2.*t27.*t32)./2.0;
t36 = (nj2.*t28.*t32)./2.0;
t37 = (nk1.*t29.*t32)./2.0;
t38 = (nj1.*t30.*t32)./2.0;
Avector = [-t35-t36-t37-t38,-nj1.*t28.*t32+nj2.*t30.*t32-nk1.*t27.*t32+nk2.*t29.*t32,vj1.*(-1.0./2.0)-vk1./2.0-(nj2.*t32.*t34)./2.0-(nk2.*t32.*t33)./2.0,t35+t36+t37+t38,vj2./2.0+vk2./2.0-(nj1.*t32.*t34)./2.0-(nk1.*t32.*t33)./2.0];
