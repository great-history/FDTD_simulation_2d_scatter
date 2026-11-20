%% 模板文件路径（原始文件）
fileName2 = 'post_processing.m';
fileName1 = 'FDTD_2D_TM_TFSF_PML_SiScatter_GPU03.m';
fileName0 = 'Dielectric_Cylinder_Setup.m';

lambdaList = [450, 500, 550, 600, 650, 700, 750]; % wavelength : nm

tic
currentDateTime = datestr(now, 'yyyymmdd_HHMMSS');

for i = 1:length(lambdaList)
    lambdaNow = lambdaList(i);

    % 读取模板文件
    txt = fileread(fileName0);
    
    % 替换 Ez 的那一行（你用正则表达式替换）
    newLine = sprintf('lambda = %.1f;', lambdaNow);
    txt = regexprep(txt, 'lambda\s*=\s*.*?;', newLine);

    % % 生成临时文件
    % tempFile = sprintf('temp_Ez_%d.m', i);

    % 直接覆盖原始文件
    fid = fopen(fileName0, 'w');
    fwrite(fid, txt);
    fclose(fid);
    
    % 运行临时文件
    run(fileName1);
    run(fileName2);
end
toc