function I = ImportRaw (filename)
    fid = fopen(filename, 'r');
    I = fread(fid, [2448, 2048], 'uint16');
    fclose(fid);
    I = I';
    I = I / 65535;
end