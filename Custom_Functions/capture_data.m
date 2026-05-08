function data = capture_data(rx, bf, bf_TDD)
    bf.Burst = false; bf.Burst = true; bf.Burst = false;
    rawdata = rx();
    data = arrangePulseData(rawdata, rx, bf, bf_TDD);
end