//======================================================================
//	Note:			Testbench for ECDH elliptic curve scalar
//	Author:			YuJia, Chen
//======================================================================

`timescale 1ns/100ps

module test_Multinv;

parameter CYCLE = 10;
parameter MAXPATITER = 512;
parameter MAXCYCLE = 600;
parameter BW_GF = 256;
parameter BW_MI = BW_GF + 1;

integer i;

reg clk;
reg en;
reg [BW_GF-1:0] a;
wire [BW_GF-1:0] value;
wire [9-1:0] power;
wire valid;

reg [10-1:0] cnt;
reg [BW_MI-1:0] data[0:6*MAXPATITER-1];

MultInv u_top(
	.clk(clk),
	.en(en),
	.a(a),
	.value(value),
	.power(power),
	.valid(valid)
);

always #(CYCLE/2) clk = ~clk;

// signal initialization
initial begin
	$display("\n\n\n>> TESTBENCH START ****************************");
	clk = 0;
	en = 0;
	$display(">> Reading data...");
	$readmemh("./pat256/multinv.dat", data);
	
	@(negedge clk);
	@(negedge clk) en = 1;
	@(negedge clk) en = 0;

end

reg err_flag;
// feeding input signals, then compare the calculation ------------------------
initial begin
	wait(en);
	
	a = 256'h62C151D9_4EE0B260_83265B9F_49837211_A840A5E2_3D96598C_89A66817_5A3AB570;
	cnt = 0;

	wait(u_top.busy);
	@(negedge clk);
	while(1) begin
		err_flag = 0;
		@(negedge clk);
		$display("round %3d", cnt);
		// if (data[cnt * 6 + 0] !== {1'b1, u_top.u}) begin
			$display("value inconsistent in round %d", cnt);
			$display("SFT u\t %65H", data[cnt * 6 + 0]);
			$display("RTL u\t %65H", u_top.u);
			// err_flag = 1;
		// end
		// if (data[cnt * 6 + 1] !== u_top.v) begin
			$display("value inconsistent in round %d", cnt);
			$display("SFT v\t %65H", data[cnt * 6 + 1]);
			$display("RTL v\t %65H", u_top.v);
			// err_flag = 1;
		// end
		// if (data[cnt * 6 + 2] !== u_top.s) begin
			$display("value inconsistent in round %d", cnt);
			$display("SFT s\t %65H", data[cnt * 6 + 2]);
			$display("RTL s\t %65H", u_top.s);
			// err_flag = 1;
		// end
		// if (data[cnt * 6 + 3] !== u_top.r) begin
			$display("value inconsistent in round %d", cnt);
			$display("SFT r\t %65H", data[cnt * 6 + 3]);
			$display("RTL r\t %65H", u_top.r);
			// err_flag = 1;
		// end
		// if (data[cnt * 6 + 4] !== u_top.x) begin
			$display("value inconsistent in round %d", cnt);
			$display("SFT x\t %65H", data[cnt * 6 + 4]);
			$display("RTL x\t %65H", u_top.x);
			// err_flag = 1;
		// end
		// if (data[cnt * 6 + 5] !== u_top.y) begin
		// 	$display("value inconsistent in round %d", cnt);
		// 	$display("SFT y\t %65H", data[cnt * 6 + 5]);
		// 	$display("RTL y\t %65H", u_top.y);
		// 	err_flag = 1;
		// end
		cnt = cnt + 1;	
		
		#5;
		if (err_flag)
			$finish;

	end
end

initial begin
	wait(valid);
	$display(">> \tu %64h", u_top.u);
	$display(">> \tv %64h", u_top.v);
	$display(">> \tr %64h", u_top.r);
	$display(">> \ts %64h", u_top.s);
	$display(">> \tk %64h", u_top.k);
	$display(">> \tx %64h", u_top.x);
	$display(">> \ty %64h", u_top.y);
	// $display(">> My Q = ")h 
	
	// if (Qx !== endpoint[0] | Qy !== endpoint[1]) begin
	// 	$display("\n>> ====================================================");
	// 	$display(">>      					ERROR");
	// 	$display(">> ====================================================\n");
	// end else begin
	// 	$display("\n>> ========!!VERIFICATION PASS, CONGRATULATION!!========\n");
	// 	$display("    ██████╗ ██████╗ ██████╗ ██████╗ ███████╗ ██████╗████████╗");
	// 	$display("   ██╔════╝██╔═══██╗██╔══██╗██╔══██╗██╔════╝██╔════╝╚══██╔══╝");
	// 	$display("   ██║     ██║   ██║██████╔╝██████╔╝█████╗  ██║        ██║   ");
	// 	$display("   ██║     ██║   ██║██╔══██╗██╔══██╗██╔══╝  ██║        ██║   ");
	// 	$display("   ╚██████╗╚██████╔╝██║  ██║██║  ██║███████╗╚██████╗   ██║   ");
	// 	$display("    ╚═════╝ ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝ ╚═════╝   ╚═╝   ");
	// 	$display("\n>> ====================================================\n");
	// end
	$finish;
end

// export waveform
initial begin
	$fsdbDumpfile("multinv.fsdb");
	$fsdbDumpvars("+mda");
end

// terminate simulation if it takes too long
initial begin
	#(CYCLE * MAXCYCLE);
	$display("\n===================================================");
	$display("      Error!!! Simulation time is too long...      ");
	$display("   There might be something wrong in your code.    ");
 	$display("===================================================\n");
 	$finish;
end

endmodule