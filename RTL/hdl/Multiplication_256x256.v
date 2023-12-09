// calculate after reset

module Multiplication_256x256 #(
    parameter BW_GF = `BW_GF
)
(
	input clk,
	input srst_n,
	input start,				// enable multiply and reset buffer
	input [BW_GF-1:0] a,	
	input [BW_GF-1:0] b,

	output [BW_GF-1:0] out,
	output valid
);

integer i;

localparam BW_SECT = `MULTLEN;
localparam COMP_CYCLE = BW_GF/BW_SECT;

reg enable, enable_n;
reg valid_prod, valid_prod_n;
reg [($clog2(COMP_CYCLE))-1:0] cnt;
reg [($clog2(COMP_CYCLE))-1:0] cnt_n;
reg [BW_SECT + BW_GF - 1:0] prod_sect;
reg [2*BW_GF-1:0] c, c_n;

reg [BW_GF-1:0] a_buf;
reg [BW_SECT-1:0] b_section[0:COMP_CYCLE-1];
reg [BW_SECT-1:0] b_buf;
reg start_buf;

always @(*)
	for (i = 0; i < COMP_CYCLE; i = i + 1)
		b_section[i] = b[i*BW_SECT +: BW_SECT];

always @(posedge clk) begin
	a_buf <= a;
	b_buf <= b_section[cnt_n];
end

always @(posedge clk)
	if (~srst_n)
		start_buf <= 0;
	else
		start_buf <= start;

always @(*)
	if (enable && cnt == COMP_CYCLE-1) begin
		valid_prod_n = 1;
	end	else begin
		valid_prod_n = 0;
	end

always @(posedge clk)
	if (~srst_n)
		valid_prod <= 0;
	else
		valid_prod <= valid_prod_n;

always @(*)
	if (start_buf) begin
		enable_n = 1;
	end else if (cnt == COMP_CYCLE-1) begin
		enable_n = 0;
	end else begin
		enable_n = enable;
	end

always @(posedge clk)
	if (~srst_n)
		enable <= 0;
	else
		enable <= enable_n;

always @* begin
	if (enable)
		cnt_n = cnt + 1;
	else 
		cnt_n = 0;
end 

always @(posedge clk)
	if (~srst_n)
		cnt <= 0;
	else
		cnt <= cnt_n;

always @(*) begin
	prod_sect = (a_buf * b_buf) + c[2*BW_GF-1 -: BW_GF];
	c_n = {prod_sect, c[BW_GF-1 -: BW_GF-BW_SECT]};
end

always @(posedge clk) begin
	if (start_buf)
		c <= 0;
	else
		c <= c_n;
end

mod_256 u_mod256_MULT(
	.clk(clk),
	.a(c),
	.valid(valid_prod),
	.b(out),
	.finish(valid)
);

endmodule
