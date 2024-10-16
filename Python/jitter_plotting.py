import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_jitter(df, genes, title=None, x_tick_fontsize=None, title_font_size = None, fig_name=None, output_path=None):
    # Group the DataFrame by Sample_ID
    #group = df.groupby('Mouse_Ear_Tag')
    group = df.groupby('Sample_ID')

    # Define a color map for distinguishing different genes
    cmap = plt.get_cmap('tab10')

    # Create subplots
    num_subplots = len(group)
    fig, axes = plt.subplots(num_subplots, 1, figsize=(17, 3*num_subplots))
    plt.suptitle(title)
        
    # Iterate through each Sample_ID and plot on its corresponding subplot
    for (sample_id, group), ax in zip(group, axes):
        # Iterate through each unique gene
        for i, gene in enumerate(genes):
            # Filter the data for the current gene
            gene_group = group[group['Numbered_gene_name'] == gene]
            
            # If the gene group is empty, plot an empty scatter point
            if gene_group.empty:
                ax.scatter(gene, 0, s=0, label=gene, color='grey') 
            else:
                # Plot the gene using sns.stripplot
                sizes = np.log10(gene_group['Cell_number'])
                sns.stripplot(x='Numbered_gene_name', y='Cell_number', data=gene_group, jitter=0.4, size=sizes, hue='Cell_number', ax=ax, legend=None)
        
        # Add labels and legend
        ax.set_ylabel('Cell Number')
        ax.set_title(sample_id)
        ax.set_xlabel(None)
        if title_font_size:
            ax.title.set_size(title_font_size)

        # Set y-axis ticks to be in 10^number
        ax.set_yscale('log')
        if x_tick_fontsize:
            ax.tick_params(axis='x', rotation=90, labelsize=x_tick_fontsize)
        else:
            ax.tick_params(axis='x', rotation=90)

    # Set x-axis label only for the last subplot
    #axes[-1].set_xlabel('Numbered_gene_name')

    # Adjust layout
    plt.tight_layout()
    
    # Save the figure if required
    if fig_name and output_path:
        fig.savefig(f'{output_path}/{fig_name}', bbox_inches='tight')



