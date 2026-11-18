// File-based testbench: reads input.hex (12-bit two's complement), writes hw_out.txt ("n y_prev y_cur").
`timescale 1ns/1ps
module tb_cfd;
    parameter integer NS  = 80000;
    parameter integer D   = 6;
    parameter integer K   = 3;
    parameter integer THR = 90;
    parameter integer HOLD = 16;

    reg clk = 0;
    reg rst = 1;
    reg signed [11:0] x = 0;
    wire valid;
    wire [31:0] n_out;
    wire signed [15:0] yp_out, yc_out;
    reg [11:0] mem [0:NS-1];
    integer i, fd;

    cfd #(.W(12), .D(D), .K(K), .THR(THR), .HOLD(HOLD)) dut (
        .clk(clk), .rst(rst), .x(x), .valid(valid),
        .n_out(n_out), .yp_out(yp_out), .yc_out(yc_out));

    always #0.5 clk = ~clk;                       // 1 GHz sampling clock

    initial begin
        $readmemh("input.hex", mem);
        fd = $fopen("hw_out.txt", "w");
        repeat (4) @(negedge clk);
        rst = 0;
        for (i = 0; i < NS; i = i + 1) begin
            x = mem[i];
            @(negedge clk);
            if (valid) $fwrite(fd, "%0d %0d %0d\n", n_out, yp_out, yc_out);
        end
        repeat (4) @(negedge clk);
        $fclose(fd);
        $finish;
    end
endmodule
