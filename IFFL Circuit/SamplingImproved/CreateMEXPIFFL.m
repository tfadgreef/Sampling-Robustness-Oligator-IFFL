home=pwd;
cd '/home/hroekel/ODEMEXv11'
parser=pwd;
cd '/home/hroekel/ODEMEXv11/Parser/CVode'
setupCVode 
cd(parser)
addPaths
for i=1:20
    eval(['mStruct.s.x',num2str(i),'=',num2str(i),';']);
    eval(['mStruct.p.par',num2str(i),'=',num2str(i),';']);
end
for i=21:49
    eval(['mStruct.p.par',num2str(i),'=',num2str(i),';']);
end
mStruct.u.u1=1;
mStruct.c.c1=1;
cd(home)
convertToC( mStruct, 'IFFL_aPcI.m' );
compileC( 'PIFFL' );