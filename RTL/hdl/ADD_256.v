// calculate when enable come in

module ADD_256 #(
    parameter BW_GF = `BW_GF
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
reg valid_mod;
reg [BW_GF:0] sub0;
reg [BW_GF:0] add0;
reg [BW_GF:0] c, c_n;
reg [BW_GF:0] c_tmp, c_tmp_n;
wire [BW_GF:0] p;

assign p = `PRIME;

always @* begin
	add0 = a + b;
	sub0 = a + (p - b);
	c_tmp_n = is_sub ? sub0 : add0;
end 

always @(posedge clk) begin
	c_tmp <= c_tmp_n;
	valid_add <= en;
end

always @(*) begin
	if (c_tmp > p)
		c_n = c_tmp - p;
	else
		c_n = c_tmp;
end

always @(posedge clk) begin
	valid_mod <= valid_add;
	c <= c_n;
end

assign out = c[`BW_GF-1:0];
assign valid = valid_mod;

endmodule
