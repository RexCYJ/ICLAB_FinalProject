// calculate after reset

module Multiplication_192x192 #(
    parameter BW_GF = `BW_GF
)
(
	input clk,
	input rst,				// enable multiply and reset buffer
	input [BW_GF-1:0] a,	
	input [BW_GF-1:0] b,

	output [BW_GF-1:0] out,
	output valid
);

integer i;

localparam BW_SECT = 16;
localparam COMP_CYCLE = BW_GF/BW_SECT;

reg enable;
reg valid_prod;
// reg [BW_SECT-1:0] b_section;
reg [($clog2(COMP_CYCLE))-1:0] cnt;
reg [($clog2(COMP_CYCLE))-1:0] cnt_n;
reg [BW_SECT + BW_GF - 1:0] prod_sect;
reg [2*BW_GF-1:0] c, c_n;

reg [BW_GF-1:0] a_buf;
reg [BW_SECT-1:0] b_buf;
reg rst_buf;

reg [BW_SECT-1:0] b_section[0:COMP_CYCLE-1];

always @(*) begin
	for (i = 0; i < COMP_CYCLE; i = i + 1)
		b_section[i] = b[i*BW_SECT +: BW_SECT];
end

always @(posedge clk) begin
	a_buf <= a;
	b_buf <= b_section[cnt_n];
	rst_buf <= rst;
end

always @(posedge clk) begin
	if (enable && cnt == COMP_CYCLE-1) begin
		valid_prod <= 1;
	end
	else begin
		valid_prod <= 0;
	end
end

always @(posedge clk) begin
	if (rst_buf) begin
		enable <= 1;
	end else if (rst | cnt == COMP_CYCLE-1) begin
		enable <= 0;
	end
end

always @* begin
	if (~enable) begin
		cnt_n = 0;
	end
	else begin
		cnt_n = cnt + 1;
	end
end 

always @(posedge clk) begin
	if (rst_buf)
		cnt <= 0;
	else
		cnt <= cnt_n;
end

always @(*) begin
	// b_section = b_buf[cnt*BW_SECT +: BW_SECT];
	prod_sect = (a_buf * b_buf) + c[2*BW_GF-1 -: BW_GF];
	c_n = {prod_sect, c[BW_GF-1 -: BW_GF-BW_SECT]};
end

always @(posedge clk) begin
	if (rst_buf)
		c <= 0;
	else
		c <= c_n;
end

mod_192 u_mod192_MULT(
	.clk(clk),
	.a(c),
	.valid(valid_prod),
	.b(out),
	.finish(valid)
);

endmodule
