//==================================================================================================
//  Note:			2022.12.31 Fall2022 NTHU ICLAB Final Project
//  Project:		AES encryption and decryption with ECC key encryption
//  Author:			Boris
//	Cycles:			2 (without output buffer)
//==================================================================================================

module mod_192 #(
	parameter a_len = 2*`BW_GF,
	parameter b_len = `BW_GF
)
(
	input clk,
	input [a_len-1:0] a,
	input valid,
	output reg [b_len-1:0] b,
	output reg finish
);

integer i;

reg [63:0] A [0:5];
wire [b_len-1:0] T;
wire [b_len-1:0] S1, S2, S3;
reg [b_len+2-1:0] temp_b, temp_b_n;
reg [b_len+2-1:0] temp_b_q1, temp_b_q2;
reg [b_len-1:0] b_n;
reg valid_tmp;

wire [`BW_GF+2-1:0] p;

assign p = {2'b0, `PRIME};

always @(*)
	for (i = 0; i < 6; i = i + 1) 
		A[i] = a[i * 64 +: 64];

assign T  = {A[2], A[1],  A[0]};
assign S1 =       {A[3],  A[3]};
assign S2 = {A[4], A[4], 64'b0};
assign S3 = {A[5], A[5],  A[5]};

always @(*) begin
	temp_b_n = T + S1 + S2 + S3;
end

always @(posedge clk) begin
	temp_b <= temp_b_n;
	valid_tmp <= valid;
end

always @(*) begin
	temp_b_q1 = temp_b - p;
	temp_b_q2 = temp_b - {p[`BW_GF+1-1:0], 1'b0};	// -2p
end

always @(*) begin
	if (temp_b < p)
		b_n = temp_b;
	else if (temp_b_q1 < p)
		b_n = temp_b_q1;
	else
		b_n = temp_b_q2;
end

always @(posedge clk) begin
	b <= b_n;
	finish <= valid_tmp;
end

endmodule
