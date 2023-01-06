// calculate when enable come in

module ADD_192 #(
    parameter BW_GF = 192
)
(
	input clk,
	input en,				// enable multiply and reset buffer
	input [BW_GF-1:0] a,	
	input [BW_GF-1:0] b,
	input is_sub,

	output [BW_GF-1:0] out,
	output valid
);

reg valid_add;
reg [BW_GF:0] sub0;
reg [BW_GF:0] add0;
reg [BW_GF:0] c, c_n;
reg [BW_GF-1:0] p;

always @* begin
	add0 = a + b;
	sub0 = a + (p - b);
	c_n = is_sub ? sub0 : add0;
end 

always @(posedge clk) begin
	valid_add <= en;
	c <= c_n;
end

always @(posedge clk) begin
	p <= `PRIME;      // 2^192-2^64-1
end

mod_192 u_mod192_ADD(
	.clk(clk),
	.a({191'd0, c}),
	.valid(valid_add),
	.b(out),
	.finish(valid)
);

/*
reg valid_tmp;

always @(posedge clk) begin
	add0 <= a + b;
	sub0 <= a + (p - b);
end 

always @(posedge clk) begin
	valid_tmp <= en;
	valid <= valid_tmp;
end

always @(posedge clk) begin
	c <= is_sub ? sub0 : add0;
end

*/

endmodule
