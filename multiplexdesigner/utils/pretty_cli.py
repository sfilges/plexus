from rich.console import Console
from rich.panel import Panel
from rich.text import Text
from rich.align import Align
from rich.columns import Columns
from rich import box

def display_welcome():
    console = Console()
    
    # Create the main title
    title = Text("🧬 MULTIPLEX PRIMER DESIGNER 🧬", style="bold cyan")
    title_panel = Panel(
        Align.center(title),
        box=box.DOUBLE,
        style="bright_blue",
        padding=(1, 2)
    )
    
    # Create feature highlights
    features = [
        "🎯 Multi-target primer design",
        "⚡ High-throughput processing", 
        "🔬 PCR optimization tools",
        "📊 Comprehensive analysis",
        "🧪 Specificity validation",
        "📈 Results visualization"
    ]
    
    feature_columns = Columns(
        [Text(feature, style="green") for feature in features],
        equal=True,
        expand=True
    )
    
    # Welcome message content
    welcome_text = Text.assemble(
        ("Welcome to the ", "white"),
        ("Multiplex Primer Development Suite", "bold yellow"),
        ("!\n\n", "white"),
        ("This powerful tool helps you design optimized primer sets for ", "white"),
        ("multiplex PCR applications", "bold green"),
        (", ensuring high specificity and efficiency across multiple targets.", "white")
    )
    
    # Instructions
    instructions = Text.assemble(
        ("🚀 ", "yellow"),
        ("Getting Started:", "bold white"),
        ("\n• Load your target sequences", "cyan"),
        ("\n• Configure primer parameters", "cyan"), 
        ("\n• Run the optimization engine", "cyan"),
        ("\n• Export your primer sets", "cyan")
    )
    
    # Display everything
    console.print()
    console.print(title_panel)
    console.print()
    
    console.print(Panel(
        welcome_text,
        title="[bold white]About[/bold white]",
        border_style="green",
        padding=(1, 2)
    ))
    
    console.print()
    console.print(Panel(
        feature_columns,
        title="[bold white]Features[/bold white]",
        border_style="yellow",
        padding=(1, 1)
    ))
    
    console.print()
    console.print(Panel(
        instructions,
        title="[bold white]Quick Start[/bold white]",
        border_style="cyan",
        padding=(1, 2)
    ))
    
    console.print()
    console.print("[dim]Press Enter to continue...[/dim]", style="italic")
    input()