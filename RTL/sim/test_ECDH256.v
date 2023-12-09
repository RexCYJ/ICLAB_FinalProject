//======================================================================
//	Note:			Testbench for ECDH elliptic curve scalar
//	Author:			YuJia, Chen
//======================================================================

`timescale 1ns/100ps

`define BASEPTR 0

module test_ECDH256;

parameter CYCLE = 10;
parameter MAXPATNUM = 100000;
parameter MAXCYCLE = 500000;
parameter BW_GF = 256;

integer i, j;
integer addstep, regindex;

integer cycle_cnt;

reg clk, srst_n;
reg start;
reg [BW_GF-1:0] k;
reg [BW_GF-1:0] Px;
reg [BW_GF-1:0] Py;
wire [BW_GF-1:0] Qx;
wire [BW_GF-1:0] Qy;
wire valid;

reg [BW_GF-1:0] basepoint[0:2*MAXPATNUM-1];
reg [BW_GF-1:0] scalar[0:MAXPATNUM-1];
reg [BW_GF-1:0] endpoint[0:2*MAXPATNUM-1];
reg [BW_GF-1:0] regfile[0:MAXPATNUM-1];
reg [BW_GF-1:0] addtrace[0:7*MAXPATNUM-1];
reg stage_result_sw;
reg pntadd_result_sw;

ECpoint_scalar u_top(
	.clk(clk),
	.rst_n(srst_n),
	.start(start),
	.k(k),
	.Px(Px),	.Py(Py),
	.Qx(Qx),	.Qy(Qy),
	.valid(valid)
);

always #(CYCLE/2) clk = ~clk;

initial begin
	`ifdef STAGE_RESULT
		stage_result_sw = 1;
	`else
		stage_result_sw = 0;
	`endif

	`ifdef PNTADD_RESULT
		pntadd_result_sw = 1;
	`else
		pntadd_result_sw = 0;
	`endif 
end


// signal initialization
initial begin
	$display("\n\n\n>> TESTBENCH START ****************************");
	clk = 0;
	srst_n = 1;
	start = 0;
	$display(">> Reading: base point");
	$readmemh("./pat256/basepoint.dat", basepoint);
	$readmemh("./pat256/endpoint.dat", endpoint);
	$readmemh("./pat256/scalar.dat", scalar);
	$readmemh("./pat256/reg_record.dat", regfile);
	$readmemh("./pat256/pointadd.dat", addtrace);
	
	@(negedge clk) srst_n = 0;
	@(negedge clk) srst_n = 1;
	@(negedge clk);

end

// feeding input signals
initial begin
	wait(~srst_n);
	wait(srst_n);
	
	@(negedge clk);
	$display("\n\n\n>> CRYPTO START ****************************");
	start = 1;
	Px = basepoint[0 + `BASEPTR*2];
	Py = basepoint[1 + `BASEPTR*2];
	k = scalar[`BASEPTR];	
	
	@(negedge clk);
	start = 0;
	cycle_cnt = 0;	
	
	while(1) begin
		@(negedge clk);
		cycle_cnt = cycle_cnt + 1;
	end
end

reg err_flag;

// investigate each end of state
initial begin
	i = 0;
	wait(start);
	while(~valid & stage_result_sw) begin
		err_flag = 0;
		@(negedge clk);
		if (u_top.done_pointcal & u_top.state !== 7) begin
			if (u_top.state == 4)
				$display("\nIn Point Doubling");
			else if (u_top.state == 5)
				$display("\nIn Point Addition");
			else if (u_top.state == 7)
				$display("\nIn MultInv");
			
			$display("Register Bank value:");

			for (j = 0; j < 4; j = j + 1) begin
				if (regfile[j + i * 4] !== u_top.data[j]) begin
					$display("!!! reg %c mismatch !!!", 65+j);
					err_flag = 1;
				end
				$display("GOLDEN: %c: %64h", 65+j, regfile[j + i * 4]);
				$display("RTL:    %c: %64h", 65+j, u_top.data[j]);
			end
			i = i + 1;
			#5;
			if (err_flag)
				$finish;
		end
	end
end


initial begin
	addstep = 0;
	regindex = 0;
	wait(~u_top.is_first);
	while(~valid & pntadd_result_sw) begin
		@(negedge clk);
		if (u_top.state == 5 & u_top.done_calc == 1 & u_top.step_cnt > 0) begin
			@(negedge clk);
			for (regindex = 0; regindex < 9; regindex = regindex + 1) begin
				if (addtrace[addstep * 9 + regindex] !== u_top.data[regindex]) begin
					$display("Point Add: %2d - %2d: reg %c is wrong!", addstep/13, addstep % 13 + 1, 65 + regindex);
					$display("Golden: %64h", addtrace[addstep * 9 + regindex]);
					$display("   RTL: %64h", u_top.data[regindex]);
					$finish;
				end
			end
			#5;
			addstep = addstep + 1;
		end
	end
end

initial begin
	wait(valid);
	$display("\n\n>> FINISH COMPUTING =========================================");
	$display(">> Golden = ");
	$display(">> \tX %64h", endpoint[0 + `BASEPTR*2]);
	$display(">> \tY %64h", endpoint[1 + `BASEPTR*2]);
	$display(">> RTL = ");
	$display(">> \tX %64h", Qx);
	$display(">> \tY %64h", Qy);
	$display(" Total cycle:  %d", cycle_cnt);
	
	if (Qx !== endpoint[0] | Qy !== endpoint[1]) begin
		$display("\n>> ===========================================================");
		$display(">>      					ERROR");
		$display(">> ===========================================================\n");
	end else begin
		$display("\n>> ========!!VERIFICATION PASS, CONGRATULATION!!===============\n");
		$display("    ██████╗ ██████╗ ██████╗ ██████╗ ███████╗ ██████╗████████╗");
		$display("   ██╔════╝██╔═══██╗██╔══██╗██╔══██╗██╔════╝██╔════╝╚══██╔══╝");
		$display("   ██║     ██║   ██║██████╔╝██████╔╝█████╗  ██║        ██║   ");
		$display("   ██║     ██║   ██║██╔══██╗██╔══██╗██╔══╝  ██║        ██║   ");
		$display("   ╚██████╗╚██████╔╝██║  ██║██║  ██║███████╗╚██████╗   ██║   ");
		$display("    ╚═════╝ ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝ ╚═════╝   ╚═╝   ");
		$display(">> ===========================================================\n");
	end
	$finish;
end

// export waveform
initial begin
	$fsdbDumpfile("ECDH256.fsdb");
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